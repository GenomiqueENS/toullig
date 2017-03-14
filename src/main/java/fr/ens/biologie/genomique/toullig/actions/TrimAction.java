package fr.ens.biologie.genomique.toullig.actions;

import java.io.File;
import java.util.List;

import fr.ens.biologie.genomique.eoulsan.actions.AbstractAction;
import fr.ens.biologie.genomique.toullig.TrimFastq;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import fr.ens.biologie.genomique.eoulsan.Common;
import fr.ens.biologie.genomique.toullig.Globals;

import static fr.ens.biologie.genomique.eoulsan.EoulsanLogger.getLogger;

/**
 * This class define the Local exec Action.
 * @since 1.0
 * @author Laurent Jourdren
 */
public class TrimAction extends AbstractAction {

    /** Name of this action. */
    public static final String ACTION_NAME = "Trim";

    //
    // Action methods
    //

    @Override
    public String getName() {
        return ACTION_NAME;
    }

    @Override
    public String getDescription() {
        return "execute " + Globals.APP_NAME + " in local mode.";
    }

    @Override
    public void action(final List<String> arguments) {

        final Options options = makeOptions();
        final CommandLineParser parser = new GnuParser();

        File samFile = new File("");
        File fastqFile=new File("");
        File fastqOutputFile=new File("");

        try {

            // parse the command line arguments
            final CommandLine line = parser.parse(options,
                    arguments.toArray(new String[arguments.size()]), true);

            // Display help
            if (line.hasOption("help")) {
                help(options);
            }

            {
                String[] remainder = line.getArgs();
                if (remainder.length >= 2) {
                    samFile = new File(remainder[0]);
                    fastqFile = new File(remainder[1]);
                    fastqOutputFile = new File(remainder[2]);
                }
            }

        } catch (ParseException e) {
            Common.errorExit(e,
                    "Error while parsing command line arguments: " + e.getMessage());
        }
        // Execute program in local mode
        run(samFile, fastqFile, fastqOutputFile);
    }

    //
    // Command line parsing
    //

    /**
     * Create options for command line
     * @return an Options object
     */
    @SuppressWarnings("static-access")
    private static Options makeOptions() {

        // create Options object
        final Options options = new Options();

        options.addOption(OptionBuilder.withArgName("help").hasArg()
                .withDescription(
                        "display help")
                .create("help"));

        options
                .addOption(OptionBuilder.withArgName("samFile").hasArg()
                        .withDescription("the path to the .sam file")
                        .create("samFile"));

        options.addOption(OptionBuilder.withArgName("fastqFile").hasArg()
                .withDescription("the path to the .fastq file")
                .create("fastqFile"));

        options.addOption(OptionBuilder.withArgName("fastqOutputFile").hasArg()
                .withDescription("the path to the .fastq file")
                .create("fastqOutputFile"));

        return options;
    }

    /**
     * Show command line help.
     * @param options Options of the software
     */
    private static void help(final Options options) {

        // Show help message
        final HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp(
                Globals.APP_NAME_LOWER_CASE
                        + ".sh " + ACTION_NAME + " [options] workflow.xml design.txt",
                options);

        Common.exit(0);
    }

    //
    // Execution
    //


    /**
     * Run Toullig Trim
     * @param samFile, a sam file
     * @param fastqFile, a fasqt file
     * @param fastqOutputFile, a fastq trimmed at output
     */
    private static void run(final File samFile, final File fastqFile, final File fastqOutputFile) {

        try {

            getLogger().info("Sam File : " + samFile);
            getLogger().info("Fastq File: " + fastqFile);
            getLogger().info("Fastq Trimmed Output File: " + fastqOutputFile);

            TrimFastq clean = new TrimFastq(samFile,fastqFile,fastqOutputFile);

        }catch (Exception e3){
            e3.printStackTrace();
        }
    }
}