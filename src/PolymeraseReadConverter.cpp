// Author: Derek Barnett

#include "PolymeraseReadConverter.h"
#include <pbbam/BamRecord.h>
#include <pbbam/BamWriter.h>

PolymeraseReadConverter::PolymeraseReadConverter(Settings& settings)
    : ConverterBase(settings)
{ }

PolymeraseReadConverter::~PolymeraseReadConverter(void) { }

bool PolymeraseReadConverter::ConvertFile(HDFBasReader* reader,
                                          PacBio::BAM::BamWriter* writer)
{
    assert(reader);

    // initialize read scores
    InitReadScores(reader);

    // fetch records from HDF5 file
    SMRTSequence smrtRecord;
    while (reader->GetNext(smrtRecord)) {

        // Skip empty records
        if ((smrtRecord.length == 0) || !IsSequencingZmw(smrtRecord))
            continue;

        // attempt convert BAX to BAM
        if (!WriteRecord(smrtRecord, 0, smrtRecord.length, ReadGroupId(), writer))
        {
            smrtRecord.Free();
            return false;
        }

        smrtRecord.Free();
    }

    // if we get here, all OK
    return true;
}

bool PolymeraseReadConverter::ConvertFile(HDFBasReader* reader,
                                          PacBio::BAM::BamWriter* writer,
                                          PacBio::BAM::BamWriter* scrapsWriter)
{ return false; }

std::string PolymeraseReadConverter::HeaderReadType(void) const
{ return "POLYMERASE"; }

std::string PolymeraseReadConverter::ScrapsReadType(void) const
{ return "UNKNOWN"; }

std::string PolymeraseReadConverter::OutputFileSuffix(void) const
{ return ".polymerase.bam"; }

std::string PolymeraseReadConverter::ScrapsFileSuffix(void) const
{ return ".empty.bam"; }
