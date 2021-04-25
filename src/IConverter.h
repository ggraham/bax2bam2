// Author: Derek Barnett

#ifndef ICONVERTER_H
#define ICONVERTER_H

#include <map>
#include <string>
#include <vector>

#include <pbbam/BamHeader.h>
#include <pbbam/BamWriter.h>
#include <pbdata/SMRTSequence.hpp>

#include "Settings.h"

namespace PacBio {
namespace BAM {

class BamRecordImpl;

} // namespace BAM
} // namespace PacBio

class IConverter
{
public:
    virtual ~IConverter(void);

public:
    virtual std::vector<std::string> Errors(void) const final;
    virtual bool Run(void) =0;

protected:
    IConverter(Settings& settings);

    virtual void AddErrorMessage(const std::string& e) final;

    virtual PacBio::BAM::BamHeader CreateHeader(const std::string& modeString) final;

    virtual std::string HeaderReadType(void) const =0;
    virtual std::string OutputFileSuffix(void) const =0;

protected:
    // common state
    Settings& settings_;
    std::vector<std::string> errors_;

    // run info for BamHeader creation
    std::string bindingKit_;
    std::string sequencingKit_;
    std::string basecallerVersion_;
    std::string frameRateHz_;
};

#endif // ICONVERTER_H
