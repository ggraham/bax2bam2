// Author: Derek Barnett

#ifndef SUBREADCONVERTER_H
#define SUBREADCONVERTER_H

#include "ConverterBase.h"

class SubreadConverter : public ConverterBase<>
{
public:
    SubreadConverter(Settings& settings);
    ~SubreadConverter(void);

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

#endif // SUBREADCONVERTER_H
