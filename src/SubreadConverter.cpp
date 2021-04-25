// Author: Derek Barnett

#include "SubreadConverter.h"

#include <algorithm>
#include <deque>
#include <memory>

#include <pbbam/BamRecord.h>
#include <pbbam/BamWriter.h>

#include <alignment/utils/RegionUtils.hpp>
#include <hdf/HDFRegionTableReader.hpp>

#define MAX( A, B )     ( (A)>(B) ? (A) : (B) )
#define MAX3( A, B, C ) MAX( MAX( A, B ), C )

using namespace PacBio;
using namespace PacBio::BAM;

SubreadConverter::SubreadConverter(Settings& settings)
    : ConverterBase(settings)
{ }

SubreadConverter::~SubreadConverter(void) { }

struct SubreadInterval
{
    size_t Start;
    size_t End;
    PacBio::BAM::LocalContextFlags LocalContextFlags;

    SubreadInterval()
        : Start{0}
        , End{0}
        , LocalContextFlags{NO_LOCAL_CONTEXT}
    { }

    SubreadInterval(size_t start, size_t end, bool adapterBefore = false, bool adapterAfter = false)
        : Start{start}
        , End{end}
        , LocalContextFlags{(adapterBefore ? ADAPTER_BEFORE : NO_LOCAL_CONTEXT) |
                            (adapterAfter  ? ADAPTER_AFTER  : NO_LOCAL_CONTEXT)}
    { }
};

namespace {

struct ReadIntervalComparer {
    bool operator()(const ReadInterval& lhs, const ReadInterval& rhs) const {
        if (lhs.start == rhs.start)
            return lhs.end < rhs.end;
        return lhs.start < rhs.start;
    }
};

SubreadInterval ComputeSubreadIntervals(std::deque<SubreadInterval>* const intervals,
                                        std::deque<SubreadInterval>* const adapters,
                                        RegionTable& regionTable,
                                        const unsigned holeNumber,
                                        const size_t readLength)
{
    constexpr int RegionStart = RegionAnnotation::REGIONSTARTCOL;
    constexpr int RegionEnd   = RegionAnnotation::REGIONENDCOL;

    // clear the input first
    intervals->clear();
    adapters->clear();

    // region annotations of a zmw
    RegionAnnotations zmwRegions = regionTable[holeNumber];

    // Has non-empty HQregion or not?
    if (!zmwRegions.HasHQRegion())
        return SubreadInterval(0, 0);

    size_t hqStart = zmwRegions.HQStart();
    size_t hqEnd   = zmwRegions.HQEnd();

    // Catch and repair 1-off errors in the HQ region
    hqEnd = (hqEnd == readLength-1) ? readLength : hqEnd;

    // Catch empty or invalid HQ regions and return empty
    if (hqEnd <= hqStart)
        return SubreadInterval(0, 0);

    // adapter intervals of this zmw
    std::vector<ReadInterval> adapterIntervals = zmwRegions.AdapterIntervals();

    // Catch and trim overlapping adapter calls
    // Shared starts indicate multiple alignments for the same adapter
    // Unique starts indicate multiple overlapping adapters
    // Therefore we trim adapter ends and remove any 0-length adapters
    // such that the number of adapter regions == number of adapters
    stable_sort(adapterIntervals.begin(), adapterIntervals.end(), ReadIntervalComparer());
    for (size_t i = 1; i < adapterIntervals.size(); i++) {
        if (adapterIntervals[i-1].end > adapterIntervals[i].start)
            adapterIntervals[i-1].end = adapterIntervals[i].start;
    }
    adapterIntervals.erase(
            std::remove_if(adapterIntervals.begin(), adapterIntervals.end(),
                    [](const ReadInterval& interval) { return interval.start == interval.end; }),
            adapterIntervals.end());

    size_t subreadStart  = hqStart;
    bool   adapterBefore = false;

    for (size_t i = 0; i < adapterIntervals.size(); i++) {

        size_t adapterStart = adapterIntervals[i].start;
        size_t adapterEnd   = adapterIntervals[i].end;

        // if we're not in the HQRegion yet, skip ahead
        if (hqStart > adapterEnd)
            continue;

        // if the adapter is beyond the HQRegion, we're done
        if (hqEnd < adapterStart)
            break;

        // If the subread is greater than length=0, save it
        if (subreadStart < adapterStart)
            intervals->emplace_back(SubreadInterval(subreadStart, adapterStart, adapterBefore, true));

        // Save the region of the adapter that overlaps the HQ region
        adapters->emplace_back(SubreadInterval(MAX3(adapterStart, hqStart, subreadStart),
                    std::min(adapterEnd, hqEnd)));

        subreadStart  = adapterEnd;
        adapterBefore = true;
    }

    // Save any region between the last adatper and the end of the HQ region as a subread
    if (subreadStart < hqEnd)
        intervals->emplace_back(SubreadInterval(subreadStart, hqEnd, adapterBefore, false));

    return SubreadInterval(hqStart, hqEnd);
}

} // anon

bool SubreadConverter::ConvertFile(HDFBasReader* reader,
                                   PacBio::BAM::BamWriter* writer)
{
    return ConvertFile(reader, writer, nullptr);
}

bool SubreadConverter::ConvertFile(HDFBasReader* reader,
                                   PacBio::BAM::BamWriter* writer,
                                   PacBio::BAM::BamWriter* scrapsWriter)
{
    assert(reader);

    // initialize with default values (shared across all unmapped subreads)
    BamRecordImpl bamRecord;

    // read region table info
    std::unique_ptr<HDFRegionTableReader> const regionTableReader(new HDFRegionTableReader);
    RegionTable regionTable;
    std::string fn = filenameForReader_[reader];
    assert(!fn.empty());
    if (regionTableReader->Initialize(fn) == 0) {
        AddErrorMessage("could not read region table on "+fn);
        return false;
    }
    regionTable.Reset();
    regionTableReader->ReadTable(regionTable);
    regionTableReader->Close();

    // initialize read scores
    InitReadScores(reader);

    // fetch records from HDF5 file
    SMRTSequence smrtRecord;
    while (reader->GetNext(smrtRecord)) {

        // compute subread & adapter intervals
        SubreadInterval hqInterval;
        std::deque<SubreadInterval> subreadIntervals;
        std::deque<SubreadInterval> adapterIntervals;
        try {
            hqInterval = ComputeSubreadIntervals(&subreadIntervals,
                                                 &adapterIntervals,
                                                 regionTable,
                                                 smrtRecord.zmwData.holeNumber,
                                                 smrtRecord.length);
        } catch (std::runtime_error& e) {
            AddErrorMessage(std::string(e.what()));
            smrtRecord.Free();
            return false;
        }

        // sequencing ZMW
        if (IsSequencingZmw(smrtRecord))
        {
            // write subreads to main BAM file
            for (const SubreadInterval& interval : subreadIntervals)
            {
                // skip invalid or 0-sized intervals
                if (interval.End <= interval.Start)
                    continue;

                if (!WriteSubreadRecord(smrtRecord,
                                        interval.Start,
                                        interval.End,
                                        ReadGroupId(),
                                        static_cast<uint8_t>(interval.LocalContextFlags),
                                        writer))
                {
                    smrtRecord.Free();
                    return false;
                }
            }

            // if scraps BAM file present
            if (scrapsWriter)
            {
                // write 5-end LQ sequence
                if (hqInterval.Start > 0)
                {
                    if (!WriteLowQualityRecord(smrtRecord,
                                               0,
                                               hqInterval.Start,
                                               ScrapsReadGroupId(),
                                               scrapsWriter))
                    {
                        smrtRecord.Free();
                        return false;
                    }
                }

                // write adapters
                for (const SubreadInterval& interval : adapterIntervals) {

                    // skip invalid or 0-sized adapters
                    if (interval.End <= interval.Start)
                        continue;

                    if (!WriteAdapterRecord(smrtRecord,
                                            interval.Start,
                                            interval.End,
                                            ScrapsReadGroupId(),
                                            scrapsWriter))
                    {
                        smrtRecord.Free();
                        return false;
                    }
                }

                // write 3'-end LQ sequence
                if (hqInterval.End < smrtRecord.length)
                {
                    if (!WriteLowQualityRecord(smrtRecord,
                                               hqInterval.End,
                                               smrtRecord.length,
                                               ScrapsReadGroupId(),
                                               scrapsWriter))
                    {
                        smrtRecord.Free();
                        return false;
                    }
                }
            }
        } // sequencing ZMW

        // non-sequencing ZMW
        else
        {
            assert(!IsSequencingZmw(smrtRecord));

            // only write these if scraps BAM present & we are in 'internal mode'
            if (settings_.isInternal && scrapsWriter)
            {
                // write 5-end LQ sequence to scraps BAM
                if (hqInterval.Start > 0)
                {
                    if (!WriteLowQualityRecord(smrtRecord,
                                               0,
                                               hqInterval.Start,
                                               ScrapsReadGroupId(),
                                               scrapsWriter))
                    {
                        smrtRecord.Free();
                        return false;
                    }
                }

                // write subreads & adapters to scraps BAM, sorted by query start
                while (!subreadIntervals.empty() && !adapterIntervals.empty()) {

                    const SubreadInterval& subread = subreadIntervals.front();
                    const SubreadInterval& adapter = adapterIntervals.front();
                    assert(subread.Start != adapter.Start);

                    if (subread.Start < adapter.Start)
                    {
                        if (!WriteFilteredRecord(smrtRecord,
                                                 subread.Start,
                                                 subread.End,
                                                 ScrapsReadGroupId(),
                                                 static_cast<uint8_t>(subread.LocalContextFlags),
                                                 scrapsWriter))
                        {
                            smrtRecord.Free();
                            return false;
                        }

                        subreadIntervals.pop_front();
                    }
                    else
                    {
                        if (!WriteAdapterRecord(smrtRecord,
                                                adapter.Start,
                                                adapter.End,
                                                ScrapsReadGroupId(),
                                                scrapsWriter))
                        {
                            smrtRecord.Free();
                            return false;
                        }
                        adapterIntervals.pop_front();
                    }
                }

                // flush any traling subread intervals
                while (!subreadIntervals.empty())
                {
                    assert(adapterIntervals.empty());
                    const SubreadInterval& subread = subreadIntervals.front();
                    if (!WriteFilteredRecord(smrtRecord,
                                             subread.Start,
                                             subread.End,
                                             ScrapsReadGroupId(),
                                             static_cast<uint8_t>(subread.LocalContextFlags),
                                             scrapsWriter))
                    {
                        smrtRecord.Free();
                        return false;
                    }

                    subreadIntervals.pop_front();
                }

                // flush any remaining adapter intervals
                while (!adapterIntervals.empty())
                {
                    assert(subreadIntervals.empty());
                    const SubreadInterval& adapter = adapterIntervals.front();
                    if (!WriteAdapterRecord(smrtRecord,
                                            adapter.Start,
                                            adapter.End,
                                            ScrapsReadGroupId(),
                                            scrapsWriter))
                    {
                        smrtRecord.Free();
                        return false;
                    }
                    adapterIntervals.pop_front();
                }

                // write 3'-end LQ sequence to scraps BAM
                if (hqInterval.End < smrtRecord.length)
                {
                    if (!WriteLowQualityRecord(smrtRecord,
                                               hqInterval.End,
                                               smrtRecord.length,
                                               ScrapsReadGroupId(),
                                               scrapsWriter))
                    {
                        smrtRecord.Free();
                        return false;
                    }
                }
            }
        } // non-sequencing ZMW

        smrtRecord.Free();
    }

    // if we get here, all OK
    return true;
}

std::string SubreadConverter::HeaderReadType(void) const
{ return "SUBREAD"; }

std::string SubreadConverter::ScrapsReadType(void) const
{ return "SCRAP"; }

std::string SubreadConverter::OutputFileSuffix(void) const
{ return ".subreads.bam"; }

std::string SubreadConverter::ScrapsFileSuffix(void) const
{ return ".scraps.bam"; }
