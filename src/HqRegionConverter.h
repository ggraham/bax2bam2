// Author: Derek Barnett

#ifndef HQREGIONCONVERTER_H
#define HQREGIONCONVERTER_H

#include "ConverterBase.h"

class HqRegionConverter : public ConverterBase<>
{
public:
    HqRegionConverter(Settings& settings);
    ~HqRegionConverter(void);

protected:
    bool ConvertFile(HDFBasReader* reader,
                     PacBio::BAM::BamWriter* writer);
    bool ConvertFile(HDFBasReader* reader,
                     PacBio::BAM::BamWriter* writer,
                     PacBio::BAM::BamWriter* scrapsWriter);
    std::string HeaderReadType(void) const;
    std::string ScrapsReadType(void) const;
    std::string OutputFileSuffix(void) const;
    std::string ScrapsFileSuffix(void) const;
};

#endif // HQREGIONCONVERTER_H
