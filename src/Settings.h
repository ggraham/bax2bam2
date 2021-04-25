// Author: Derek Barnett
#ifndef SETTINGS_H
#define SETTINGS_H

#include <string>
#include <vector>

namespace optparse { class OptionParser; }

class Settings
{
public:
    enum Mode { SubreadMode
              , HQRegionMode
              , PolymeraseMode
              , CCSMode
              };

    struct Option {
        static const char* datasetXml_;
        static const char* hqRegionMode_;
        static const char* input_;
        static const char* fofn_;
        static const char* losslessFrames_;
        static const char* output_;
        static const char* polymeraseMode_;
        static const char* pulseFeatures_;
        static const char* subreadMode_;
        static const char* ccsMode_;
        static const char* internalMode_;
        static const char* outputXml_;
        static const char* sequelPlatform_;
        static const char* allowUnsupportedChem_;
    };

public:
    Settings(void);
    static Settings FromCommandLine(optparse::OptionParser& parser,
                                    int argc,
                                    char* argv[]);

public:
    // input/output
    std::vector<std::string> inputFilenames;
    std::vector<std::string> inputBaxFilenames;
    std::string datasetXmlFilename;
    std::string fofnFilename;
    std::string outputBamPrefix;
    std::string outputBamFilename;
    std::string scrapsBamFilename;
    std::string outputXmlFilename;

    // mode
    Mode mode;
    bool isInternal;

    // platform
    bool isSequelInput;

    // chemistry checking?
    bool isIgnoringChemistryCheck;

    // features
    bool usingDeletionQV;
    bool usingDeletionTag;
    bool usingInsertionQV;
    bool usingIPD;
    bool usingMergeQV;
    bool usingPulseWidth;
    bool usingSubstitutionQV;
    bool usingSubstitutionTag;

    // frame data encoding
    bool losslessFrames;

    // program info
    std::string program;
    std::string args;
    std::string version;
    std::string description;

    // generated
    std::string movieName;
    std::string readGroupId;
    std::string scrapsReadGroupId;

    // command line parsing
    std::vector<std::string> errors;
};

#endif // SETTINGS_H
