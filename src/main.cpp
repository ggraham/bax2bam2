// Author: Derek Barnett

#include "Bax2Bam.h"
#include "OptionParser.h"
#include "Settings.h"
#include <iostream>
#include <string>
#include <cstdlib>

int main(int argc, char* argv[])
{
    // setup help & options
    optparse::OptionParser parser;
    parser.description("bax2bam converts the legacy PacBio basecall format (bax.h5) into the BAM basecall format.");
    parser.prog("bax2bam");
    parser.version("0.0.11");
    parser.add_version_option(true);
    parser.add_help_option(true);

    auto ioGroup = optparse::OptionGroup(parser, "Input/output files");
    ioGroup.add_option("")
           .dest(Settings::Option::input_)
	   .metavar("movie.1.bax.h5 movie.2.bax.h5 ...")
           .help("Input files which should be from the same movie");
    ioGroup.add_option("--xml")
           .dest(Settings::Option::datasetXml_)
           .metavar("STRING")
           .help("DataSet XML file containing a list of movie names");
    ioGroup.add_option("-f", "--fofn")
           .dest(Settings::Option::fofn_)
           .metavar("STRING")
           .help("File-of-file-names containing a list of input files");
    ioGroup.add_option("-o")
           .dest(Settings::Option::output_)
	   .metavar("STRING")
           .help("Prefix of output filenames. Movie name will be used if no prefix provided");
    ioGroup.add_option("--output-xml")
           .dest(Settings::Option::outputXml_)
           .metavar("STRING")
           .help("Explicit output XML name. If none provided via this arg, bax2bam will use -o prefix (<prefix>.dataset.xml). "
                 "If that is not specified either, the output XML filename will be <moviename>.dataset.xml");
    parser.add_option_group(ioGroup);

    auto platformGroup = optparse::OptionGroup(parser, "Input sequencing platform");
    platformGroup.add_option("--sequel-input")
                 .dest(Settings::Option::sequelPlatform_)
                 .action("store_true")
                 .help("Specify that input data is from Sequel. "
                       "bax2bam will assume RS unless this option is specified");

    auto readModeGroup = optparse::OptionGroup(parser, "Output read types (mutually exclusive)");
    readModeGroup.add_option("--subread")
                 .dest(Settings::Option::subreadMode_)
                 .action("store_true")
                 .help("Output subreads (default)");
    readModeGroup.add_option("--hqregion")
                 .dest(Settings::Option::hqRegionMode_)
                 .action("store_true")
                 .help("Output HQ regions");
    readModeGroup.add_option("--polymeraseread")
                 .dest(Settings::Option::polymeraseMode_)
                 .action("store_true")
                 .help("Output full polymerase read");
    readModeGroup.add_option("--ccs")
                 .dest(Settings::Option::ccsMode_)
                 .action("store_true")
                 .help("Output CCS sequences (requires ccs.h5 input)");
    parser.add_option_group(readModeGroup);

    auto featureGroup = optparse::OptionGroup(parser, "Pulse feature options");
    featureGroup.group_description("Configure pulse features in the output BAM. Supported features include:\n"
                                   "    Pulse Feature:    BAM tag:  Default:\n"
                                   "    DeletionQV        dq        Y\n"
                                   "    DeletionTag       dt        Y\n"
                                   "    InsertionQV       iq        Y\n"
                                   "    IPD               ip        Y\n"
                                   "    PulseWidth        pw        Y*\n"
                                   "    MergeQV           mq        Y\n"
                                   "    SubstitutionQV    sq        Y\n"
                                   "    SubstitutionTag   st        N\n"
                                   "\n"
                                   "* - PulseWidth will be ignored in CCS mode.\n"
                                   "\n"
                                   "If this option is used, then only those features listed will be included, "
                                   "regardless of the default state."
                                   );
    featureGroup.add_option("--pulsefeatures")
                .dest(Settings::Option::pulseFeatures_)
                .metavar("STRING")
                .help("Comma-separated list of desired pulse features, using the names listed above.\n");
    featureGroup.add_option("--losslessframes")
                .dest(Settings::Option::losslessFrames_)
                .action("store_true")
                .help("Store full, 16-bit IPD/PulseWidth data, instead of (default) downsampled, 8-bit encoding.");
    parser.add_option_group(featureGroup);

    auto bamModeGroup = optparse::OptionGroup(parser, "Output BAM file type");
    bamModeGroup.add_option("--internal")
                .dest(Settings::Option::internalMode_)
                .action("store_true")
                .help("Output BAMs in internal mode. Currently this indicates that "
                      "non-sequencing ZMWs should be included in the output scraps "
                      "BAM file, if applicable."
                      );
    parser.add_option_group(bamModeGroup);

    auto additionalGroup = optparse::OptionGroup(parser, "Additional options");
    additionalGroup.add_option("--allowUnrecognizedChemistryTriple")
                   .dest(Settings::Option::allowUnsupportedChem_)
                   .action("store_true")
                   .help("By default, bax2bam only allows the conversion of files "
                         "with chemistries that are supported in SMRT Analysis 3. "
                         "Set this flag to disable the strict check and allow "
                         "generation of BAM files containing legacy chemistries.");
    parser.add_option_group(additionalGroup);

    // parse command line
    Settings settings = Settings::FromCommandLine(parser, argc, argv);
    if (!settings.errors.empty()) {
        std::cerr << std::endl;
        for (const auto e : settings.errors)
            std::cerr << "ERROR: " << e << std::endl;
        std::cerr << std::endl;
        parser.print_help();
        return EXIT_FAILURE;
    }

    // main conversion
    return Bax2Bam::Run(settings);
}
