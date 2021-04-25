// Author: Derek Barnett

#include "Settings.h"
#include "OptionParser.h"

#include <sstream>

#include <boost/algorithm/string.hpp>

#include <pbbam/DataSet.h>

#include <hdf/HDFNewBasReader.hpp>

namespace internal {

static
std::vector<std::string> BaxFilenamesFromXml(const std::string& xmlFilename)
{
    using namespace PacBio::BAM;

    try {
        std::vector<std::string> filenames;

        DataSet dataset(xmlFilename);
        const std::vector<std::string> resources = dataset.ResolvedResourceIds();
        for (const std::string& resource : resources) {
            std::cerr << resource << std::endl;
            const boost::iterator_range<std::string::const_iterator> baxFound = boost::algorithm::ifind_first(resource, ".bax.h5");
            if (!baxFound.empty())
                filenames.push_back(resource);
        }
        return filenames;

    } catch (std::exception&) {
        // TODO: report error
        return std::vector<std::string>();
    }
}

static
std::vector<std::string> FilenamesFromFofn(const std::string& fileName)
{
    std::vector<std::string> retval;
    std::ifstream in_stream;
    std::string line;

    in_stream.open(fileName);

    while(!in_stream.eof())
    {
        in_stream >> line;
        if (!line.empty())
            retval.push_back(line);
        line.clear();
    }

    return retval;
}

static
bool isBasH5(const std::string& fileName)
{
    return boost::ends_with(boost::to_lower_copy(fileName), ".bas.h5");
}

static
void H5FilenamesFromBasH5(const std::string& basFileName,
                          std::vector<std::string>* const output)
{
    HDFNewBasReader reader;
    if (reader.Initialize(basFileName))
        for (const auto& baxFileName : reader.GetBaxFileNames())
            output->push_back(baxFileName);
    else
        output->push_back(basFileName);
}

} // namespace internal

// option names
const char* Settings::Option::datasetXml_     = "datasetXml";
const char* Settings::Option::hqRegionMode_   = "hqRegionMode";
const char* Settings::Option::input_          = "input";
const char* Settings::Option::fofn_           = "fofn";
const char* Settings::Option::losslessFrames_ = "losslessFrames";
const char* Settings::Option::output_         = "output";
const char* Settings::Option::polymeraseMode_ = "polymeraseMode";
const char* Settings::Option::pulseFeatures_  = "pulseFeatures";
const char* Settings::Option::subreadMode_    = "subreadMode";
const char* Settings::Option::ccsMode_        = "ccsMode";
const char* Settings::Option::internalMode_   = "internalMode";
const char* Settings::Option::outputXml_      = "outputXml";
const char* Settings::Option::sequelPlatform_ = "sequelPlatform";
const char* Settings::Option::allowUnsupportedChem_  = "allowUnsupportedChem";

Settings::Settings(void)
    : mode(Settings::SubreadMode)
    , isInternal(false)
    , isSequelInput(false)
    , isIgnoringChemistryCheck(false)
    , usingDeletionQV(true)
    , usingDeletionTag(true)
    , usingInsertionQV(true)
    , usingIPD(true)
    , usingMergeQV(true)
    , usingPulseWidth(true)
    , usingSubstitutionQV(true)
    , usingSubstitutionTag(false)
    , losslessFrames(false)
{ }

Settings Settings::FromCommandLine(optparse::OptionParser& parser,
                                   int argc,
                                   char *argv[])
{
    Settings settings;

    // general program info
    settings.program = parser.prog();
    settings.description = parser.description();
    settings.version = parser.version();
    for (int i = 1; i < argc; ++i) {
        settings.args.append(argv[i]);
        settings.args.append(" ");
    }

    const optparse::Values options = parser.parse_args(argc, argv);

    // output prefix
    // TODO: output dir ??
    settings.outputBamPrefix = options[Settings::Option::output_];
    settings.outputXmlFilename = options[Settings::Option::outputXml_];

    // input files from dataset XML ?
    if ( options.is_set(Settings::Option::datasetXml_) ) {
        settings.datasetXmlFilename = options[Settings::Option::datasetXml_];
        settings.inputBaxFilenames = internal::BaxFilenamesFromXml(settings.datasetXmlFilename);
    }

    // input files from fofn ?
    else if ( options.is_set(Settings::Option::fofn_))
    {
        settings.fofnFilename = options[Settings::Option::fofn_];
        settings.inputFilenames = internal::FilenamesFromFofn(settings.fofnFilename);
    }

    // else input files command-line args
    else
        settings.inputFilenames = parser.args();

    // Process input files to convert Bas.H5 --> Bax.h5 as needed
    for (const std::string& fn : settings.inputFilenames)
    {
        if (internal::isBasH5(fn))
            internal::H5FilenamesFromBasH5(fn, &settings.inputBaxFilenames);
        else
            settings.inputBaxFilenames.push_back(fn);
    }

    if (settings.inputBaxFilenames.empty())
        settings.errors.push_back("missing input BAX files.");

    // mode
    const bool isSubreadMode =
            options.is_set(Settings::Option::subreadMode_) ? options.get(Settings::Option::subreadMode_)
                                                           : false;
    const bool isHQRegionMode =
            options.is_set(Settings::Option::hqRegionMode_) ? options.get(Settings::Option::hqRegionMode_)
                                                            : false;
    const bool isPolymeraseMode =
            options.is_set(Settings::Option::polymeraseMode_) ? options.get(Settings::Option::polymeraseMode_)
                                                              : false;
    const bool isCCS =
            options.is_set(Settings::Option::ccsMode_) ? options.get(Settings::Option::ccsMode_)
                                                       : false;

    int modeCount = 0;
    if (isSubreadMode)    ++modeCount;
    if (isHQRegionMode)   ++modeCount;
    if (isPolymeraseMode) ++modeCount;
    if (isCCS)            ++modeCount;

    if (modeCount == 0)
        settings.mode = Settings::SubreadMode;
    else if (modeCount == 1) {
        if (isSubreadMode)    settings.mode = Settings::SubreadMode;
        if (isHQRegionMode)   settings.mode = Settings::HQRegionMode;
        if (isPolymeraseMode) settings.mode = Settings::PolymeraseMode;
        if (isCCS)            settings.mode = Settings::CCSMode;
    }
    else
        settings.errors.push_back("multiple modes selected");

    // internal file mode
    settings.isInternal = options.is_set(Settings::Option::internalMode_) ? options.get(Settings::Option::internalMode_)
                                                                          : false;

    // strict/relaxed chemistry check
    settings.isIgnoringChemistryCheck = options.is_set(Settings::Option::allowUnsupportedChem_) ? options.get(Settings::Option::allowUnsupportedChem_)
                                                                                                : false;

    // platform
    settings.isSequelInput = options.is_set(Settings::Option::sequelPlatform_) ? options.get(Settings::Option::sequelPlatform_)
                                                                                : false;

    // frame data encoding
    settings.losslessFrames = options.is_set(Settings::Option::losslessFrames_) ? options.get(Settings::Option::losslessFrames_)
                                                                                : false;

    // pulse features list
    if (options.is_set(Settings::Option::pulseFeatures_)) {

        // ignore defaults
        settings.usingDeletionQV = false;
        settings.usingDeletionTag = false;
        settings.usingInsertionQV = false;
        settings.usingIPD = false;
        settings.usingMergeQV = false;
        settings.usingPulseWidth = false;
        settings.usingSubstitutionQV = false;
        settings.usingSubstitutionTag = false;

        // apply user-requested features
        std::stringstream stream(options[Settings::Option::pulseFeatures_]);
        std::string feature;
        while(std::getline(stream, feature, ',')) {
            if      (feature == "DeletionQV")      settings.usingDeletionQV = true;
            else if (feature == "DeletionTag")     settings.usingDeletionTag = true;
            else if (feature == "InsertionQV")     settings.usingInsertionQV = true;
            else if (feature == "IPD")             settings.usingIPD = true;
            else if (feature == "MergeQV")         settings.usingMergeQV = true;
            else if (feature == "PulseWidth")      settings.usingPulseWidth = true;
            else if (feature == "SubstitutionQV")  settings.usingSubstitutionQV = true;
            else if (feature == "SubstitutionTag") settings.usingSubstitutionTag = true;
            else
                settings.errors.push_back(std::string("unknown pulse feature: ") + feature);
        }
    }

    // always disable PulseWidth tag in CCS mode
    if (isCCS)
        settings.usingPulseWidth = false;

#ifdef DEBUG_SETTINGS

    std::string modeString;
    if (settings.mode == Settings::SubreadMode)
        modeString = "subread";
    else if (settings.mode == Settings::HQRegionMode)
        modeString = "hqRegion";
    else if (settings.mode == Settings::PolymeraseMode)
        modeString = "polymerase";
    else
        modeString = "ccs";

    std::string platformString = settings.isSequelInput_ ? "Sequel" : "RS";

    std::cerr << "CommandLine: " << settings.program << " " << settings.args << std::endl
         << "Description: " << settings.description << std::endl
         << "Version:     " << settings.version << std::endl
         << "Mode:        " << modeString << std::endl
         << "Platform:    " << platformString << std::endl
         << "DeletionQV?:      " << ( settings.usingDeletionQV ? "yes" : "no" ) << std::endl
         << "DeletionTag?:     " << ( settings.usingDeletionTag ? "yes" : "no" ) << std::endl
         << "InsertionQV?:     " << ( settings.usingInsertionQV ? "yes" : "no" ) << std::endl
         << "IPD?:             " << ( settings.usingMergeQV ? "yes" : "no" ) << std::endl
         << "MergeQV?:         " << ( settings.usingIPD ? "yes" : "no" ) << std::endl
         << "PulseWidth?:      " << ( settings.usingPulseWidth ? "yes" : "no" ) << std::endl
         << "SubstitutionQV?:  " << ( settings.usingSubstitutionQV ? "yes" : "no" ) << std::endl
         << "SubstitutionTag?: " << ( settings.usingSubstitutionTag ? "yes" : "no" ) << std::endl;
#endif

    return settings;
}
