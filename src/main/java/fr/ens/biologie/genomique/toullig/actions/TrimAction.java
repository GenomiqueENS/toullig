package fr.ens.biologie.genomique.toullig.actions;

import java.io.File;
import java.util.List;

import com.google.common.primitives.Ints;
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
 * @author Aurélien Birer
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

        String trimmer="cutadapt";
        String mode="P";
        String stats = "false";
        int minlen = 100;
        double errorRateCutadapt=0.5;
        double thresholdSW=0.8;
        int lengthWindowsSW=15;
        int seedMismatchesTrimmomatic=17;
        int palindromeClipThresholdTrimmomatic =30;
        int simpleClipThreshold=7;

        File samFile = new File("");
        File fastqFile=new File("");
        File fastqOutputFile=new File("");
        File adaptorFile = new File("");
        File workDir = new File("");

        try {

            // parse the command line arguments
            final CommandLine line = parser.parse(options,
                    arguments.toArray(new String[arguments.size()]), true);

            // Display help
            if (line.hasOption("help")) {
                help(options);
            }

            if (line.hasOption("trimmer")) {
                trimmer = line.getOptionValue("trimmer").toLowerCase();
            }

            if (line.hasOption("mode")) {
                mode = line.getOptionValue("mode").toLowerCase();
            }

            if (line.hasOption("stats")) {
                stats = line.getOptionValue("stats").toLowerCase();
            }

            if (line.hasOption("minlen")) {
                minlen = Integer.parseInt(line.getOptionValue("minlen").toLowerCase());
            }

            if (line.hasOption("errorRateCutadapt")) {
                errorRateCutadapt = Long.parseLong(line.getOptionValue("errorRateCutadapt").toLowerCase());
            }

            if (line.hasOption("thresholdSW")) {
                thresholdSW = Long.parseLong(line.getOptionValue("thresholdSW").toLowerCase());
            }

            if (line.hasOption("lengthWindowsSW")) {
                lengthWindowsSW = Integer.parseInt(line.getOptionValue("lengthWindowsSW").toLowerCase());
            }

            if (line.hasOption("seedMismatchesTrimmomatic")) {
                seedMismatchesTrimmomatic = Integer.parseInt(line.getOptionValue("seedMismatchesTrimmomatic").toLowerCase());
            }

            if (line.hasOption("palindromeClipThresholdTrimmomatic")) {
                palindromeClipThresholdTrimmomatic = Integer.parseInt(line.getOptionValue("palindromeClipThresholdTrimmomatic").toLowerCase());
            }

            if (line.hasOption("simpleClipThreshold")) {
                simpleClipThreshold = Integer.parseInt(line.getOptionValue("simpleClipThreshold").toLowerCase());
            }

            {
                String[] remainder = line.getArgs();
                if (remainder.length >= 4) {
                    samFile = new File(remainder[0]);
                    fastqFile = new File(remainder[1]);
                    fastqOutputFile = new File(remainder[2]);
                    adaptorFile = new File(remainder[3]);
                    workDir = new File(remainder[4]);
                }
            }

        } catch (ParseException e) {
            Common.errorExit(e,
                    "Error while parsing command line arguments: " + e.getMessage());
        }
        // Execute program in local mode
        run(trimmer, mode, stats, minlen, errorRateCutadapt, thresholdSW, lengthWindowsSW, seedMismatchesTrimmomatic, palindromeClipThresholdTrimmomatic,  simpleClipThreshold, samFile, fastqFile, fastqOutputFile, adaptorFile, workDir);
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
                .addOption(OptionBuilder.withArgName("trimmer").hasArg()
                        .withDescription("name of trimmer use [cutadapt | trimmomatic] (default : cutadapt)")
                        .create("trimmer"));

        options
                .addOption(OptionBuilder.withArgName("mode").hasArg()
                        .withDescription("mode of trimming use [P | SW] (default : P)")
                        .create("mode"));
        options
                .addOption(OptionBuilder.withArgName("stats").hasArg()
                        .withDescription("make somes stats on the trimming [true | false] (default : false)")
                        .create("stats"));

        options
                .addOption(OptionBuilder.withArgName("minlen").hasArg()
                        .withDescription("minimum length of trimmed fastq write (default : 100)")
                        .create("minlen"));
        options
                .addOption(OptionBuilder.withArgName("errorRateCutadapt").hasArg()
                        .withDescription("error rate for cutadapt (default: 0.5")
                        .create("errorRateCutadapt"));
        options
                .addOption(OptionBuilder.withArgName("thresholdSW").hasArg()
                        .withDescription("threshold for Side-Windows processus (default: 0.8)")
                        .create("thresholdSW"));
        options
                .addOption(OptionBuilder.withArgName("seedMismatchesTrimmomatic").hasArg()
                        .withDescription("seed mismatches option for Trimmomatic (default: 17)")
                        .create("seedMismatchesTrimmomatic"));
        options
                .addOption(OptionBuilder.withArgName("palindromeClipThresholdTrimmomatic").hasArg()
                        .withDescription("palindrome clip threshold option for Trimmomatic (default: 30)")
                        .create("palindromeClipThresholdTrimmomatic"));
        options
                .addOption(OptionBuilder.withArgName("simpleClipThreshold").hasArg()
                        .withDescription("simple clip threshold option for Trimmomatic (default : 7)")
                        .create("simpleClipThreshold"));

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

        options.addOption(OptionBuilder.withArgName("adaptorFile").hasArg()
                .withDescription("the path to the adaptor file")
                .create("adaptorFile"));

        options.addOption(OptionBuilder.withArgName("workDir").hasArg()
                .withDescription("the work directory")
                .create("workDir"));


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
    private static void run(final String trimmer, final String mode, final String stats, final int minlen, final double errorRateCutadapt, final double thresholdSW, final int lengthWindowsSW, final int seedMismatchesTrimmomatic, final int palindromeClipThresholdTrimmomatic, final int simpleClipThreshold, final File samFile, final File fastqFile, final File fastqOutputFile, final File adaptorFile, final File workDir) {

        try {

            getLogger().info("Sam File : " + samFile);
            getLogger().info("Fastq File: " + fastqFile);
            getLogger().info("Fastq Trimmed Output File: " + fastqOutputFile);
            getLogger().info("Adaptor File: " + adaptorFile);
            getLogger().info("Work Directory: " + workDir);

            TrimFastq trim = new TrimFastq(samFile,fastqFile, adaptorFile, fastqOutputFile, workDir);

            if (trimmer.contains("cutadapt")) {
                trim.setProcessCutadapt(true);
            }
            if (trimmer.contains("trimmomatic")) {
                trim.setProcessTrimmomatic(true);
            }
            if (mode.contains("P")) {
                trim.setProcessPTrim(true);
            }
            if (mode.contains("SW")) {
                trim.setProcessSWTrim(true);
            }
            if (stats.contains("true")) {
                trim.setProcessStats(true);
            }

            trim.setThresholdSW(thresholdSW);

            trim.setLengthWindowsSW(lengthWindowsSW);

            trim.setErrorRateCutadapt(errorRateCutadapt);

            trim.setSeedMismatchesTrimmomatic(seedMismatchesTrimmomatic);

            trim.setPalindromeClipThresholdTrimmomatic(palindromeClipThresholdTrimmomatic);

            trim.setSimpleClipThreshold(simpleClipThreshold);

            trim.setMinimunLengthToWrite(minlen);

            trim.execution();


        }catch (Exception e3){
            e3.printStackTrace();
        }
    }
}