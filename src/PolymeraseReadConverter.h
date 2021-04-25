// Author: Derek Barnett

#ifndef POLYMERASEREADCONVERTER_H
#define POLYMERASEREADCONVERTER_H

#include "ConverterBase.h"

class PolymeraseReadConverter : public ConverterBase<>
{
public:
    PolymeraseReadConverter(Settings& settings);
    ~PolymeraseReadConverter(void);

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

#endif // POLYMERASEREADCONVERTER_H
