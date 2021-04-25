// Author: Derek Barnett

#include "CcsConverter.h"

#include <algorithm>
#include <iostream>

#include <pbbam/BamRecord.h>
#include <pbbam/BamWriter.h>

#include <alignment/utils/RegionUtils.hpp>
#include <hdf/HDFRegionTableReader.hpp>

using namespace PacBio;
using namespace PacBio::BAM;

CcsConverter::CcsConverter(Settings& settings)
    : ConverterBase(settings)
{
    settings_.usingMergeQV         = false;
    settings_.usingDeletionTag     = false;
    settings_.usingSubstitutionTag = false;
    settings_.usingIPD             = false;
    settings_.usingPulseWidth      = false;
}

CcsConverter::~CcsConverter(void) { }

bool CcsConverter::ConvertFile(HdfCcsReader* reader,
                               BamWriter* writer)
{
    assert(reader);

    // initialize read scores
    InitReadScores(reader);

    // fetch records from HDF5 file
    CCSSequence smrtRecord;
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

bool CcsConverter::ConvertFile(HdfCcsReader* reader,
                               PacBio::BAM::BamWriter* writer,
                               PacBio::BAM::BamWriter* scrapsWriter)
{ return false; }

void CcsConverter::SetSequenceAndQualities(PacBio::BAM::BamRecordImpl* bamRecord,
                                           const CCSSequence& smrtRead,
                                           const int start,
                                           const int length)
{
    recordSequence_.assign((const char*)smrtRead.seq + start, length);
    if (smrtRead.qual.Empty())
        bamRecord->SetSequenceAndQualities(recordSequence_);
    else
    {
        recordQVs_.assign((uint8_t*)smrtRead.qual.data + start,
                          (uint8_t*)smrtRead.qual.data + start + length);
        bamRecord->SetSequenceAndQualities(recordSequence_, recordQVs_.Fastq());
    }
}

void CcsConverter::AddRecordName(PacBio::BAM::BamRecordImpl* bamRecord,
                                 const UInt holeNumber,
                                 const int start,
                                 const int end)
{
    const std::string name = settings_.movieName + "/"
                      + std::to_string(holeNumber) + "/ccs";
    bamRecord->Name(name);
}

void CcsConverter::AddModeTags(PacBio::BAM::TagCollection* tags,
                               const CCSSequence& smrtRead,
                               const int start,
                               const int end)
{
    (*tags)["np"] = static_cast<int32_t>(smrtRead.numPasses);
}

CcsConverter::HdfCcsReader* CcsConverter::InitHdfReader()
{
    HdfCcsReader* reader = ConverterBase<CCSSequence, HdfCcsReader>::InitHdfReader();
    // set the reader to CCS mode
    reader->SetReadBasesFromCCS();
    return reader;
}

std::string CcsConverter::HeaderReadType(void) const
{ return "CCS"; }

std::string CcsConverter::ScrapsReadType(void) const
{ return "UNKNOWN"; }

std::string CcsConverter::OutputFileSuffix(void) const
{ return ".ccs.bam"; }

std::string CcsConverter::ScrapsFileSuffix(void) const
{ return ".empty.bam"; }
