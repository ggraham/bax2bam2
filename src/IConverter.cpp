// Author: Derek Barnett

#include "IConverter.h"
#include <pbbam/BamRecord.h>
#include <algorithm>
#include <iostream>
#include <set>
#include <cassert>
#include <cmath>

using namespace PacBio;
using namespace PacBio::BAM;

IConverter::IConverter(Settings& settings)
    : settings_(settings)
{ }

IConverter::~IConverter(void) { }

void IConverter::AddErrorMessage(const std::string& e)
{ errors_.push_back(e); }

BamHeader IConverter::CreateHeader(const std::string& modeString)
{
    BamHeader header;

    // @RG ID: <read group ID>
    //     DS: READTYPE=<HQREGION|POLYMERASE|SUBREAD>[;<Tag Manifest>;BINDINGKIT=<foo>;SEQUENCINGKIT=<bar>;BASECALLERVERSION=<42>]
    //     PL: PACBIO
    //     PU: <movieName>
    //
    const PlatformModelType platform = settings_.isSequelInput ? PlatformModelType::SEQUEL
                                                                : PlatformModelType::RS;
    ReadGroupInfo rg(settings_.movieName, modeString, platform);
    rg.BindingKit(bindingKit_)
      .SequencingKit(sequencingKit_)
      .BasecallerVersion(basecallerVersion_)
      .FrameRateHz(frameRateHz_);

    if (settings_.usingDeletionQV)      rg.BaseFeatureTag(BaseFeature::DELETION_QV,      "dq");
    if (settings_.usingDeletionTag)     rg.BaseFeatureTag(BaseFeature::DELETION_TAG,     "dt");
    if (settings_.usingInsertionQV)     rg.BaseFeatureTag(BaseFeature::INSERTION_QV,     "iq");
    if (settings_.usingMergeQV)         rg.BaseFeatureTag(BaseFeature::MERGE_QV,         "mq");
    if (settings_.usingSubstitutionQV)  rg.BaseFeatureTag(BaseFeature::SUBSTITUTION_QV,  "sq");
    if (settings_.usingSubstitutionTag) rg.BaseFeatureTag(BaseFeature::SUBSTITUTION_TAG, "st");
    if (settings_.usingIPD) {
        FrameCodec codec = FrameCodec::V1;
        if (settings_.losslessFrames)
            codec = FrameCodec::RAW;
        rg.IpdCodec(codec, "ip");
    }
    if (settings_.usingPulseWidth) {
        FrameCodec codec = FrameCodec::V1;
        if (settings_.losslessFrames)
            codec = FrameCodec::RAW;
        rg.PulseWidthCodec(codec, "pw");
    }

    header.AddReadGroup(rg);

    // @PG ID:bax2bam-<version>
    //     PN:bax2bam
    //     CL:bax2bam <args>
    //     DS:<description>
    //     VN:<version>

    ProgramInfo program(settings_.program + "-" + settings_.version);
    program.Name(settings_.program)
           .CommandLine(settings_.program + " " + settings_.args)
           .Description(settings_.description)
           .Version(settings_.version);
    header.AddProgram(program);

    return header;
}

std::vector<std::string> IConverter::Errors(void) const
{ return errors_; }
