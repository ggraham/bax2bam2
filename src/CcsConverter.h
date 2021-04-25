// Author: Derek Barnett

#ifndef CCSCONVERTER_H
#define CCSCONVERTER_H

#include <hdf/HDFCCSReader.hpp>
#include <pbdata/CCSSequence.hpp>

#include "ConverterBase.h"

class CcsConverter : public ConverterBase<CCSSequence, HDFCCSReader<CCSSequence>>
{
private:
    typedef HDFCCSReader<CCSSequence> HdfCcsReader;

public:
    CcsConverter(Settings& settings);
    ~CcsConverter(void);

protected:
    bool ConvertFile(HdfCcsReader* reader,
                     PacBio::BAM::BamWriter* writer);
    bool ConvertFile(HdfCcsReader* reader,
                     PacBio::BAM::BamWriter* writer,
                     PacBio::BAM::BamWriter* scrapsWriter);
    void SetSequenceAndQualities(PacBio::BAM::BamRecordImpl* bamRecord,
                                 const CCSSequence& smrtRecord,
                                 const int start,
                                 const int end);
    void AddRecordName(PacBio::BAM::BamRecordImpl* bamRecord,
                       const UInt holeNumber,
                       const int start,
                       const int end);
    void AddModeTags(PacBio::BAM::TagCollection* tags,
                     const CCSSequence& smrtRecord,
                     const int start,
                     const int end);
    HdfCcsReader* InitHdfReader(void);
    std::string HeaderReadType(void) const;
    std::string ScrapsReadType(void) const;
    std::string OutputFileSuffix(void) const;
    std::string ScrapsFileSuffix(void) const;

protected:
    PacBio::BAM::QualityValues recordQVs_;
};

#endif // CCSCONVERTER_H
