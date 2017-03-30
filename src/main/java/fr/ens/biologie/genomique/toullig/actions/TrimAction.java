package fr.ens.biologie.genomique.toullig.actions;

import static fr.ens.biologie.genomique.eoulsan.EoulsanLogger.getLogger;

import java.io.File;
import java.util.List;

import org.apache.commons.cli.*;

import fr.ens.biologie.genomique.eoulsan.Common;
import fr.ens.biologie.genomique.eoulsan.actions.AbstractAction;
import fr.ens.biologie.genomique.toullig.Globals;
import fr.ens.biologie.genomique.toullig.TrimFastq;

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

    String trimmer = "";
    String mode = "";
    String stats = "";
    double errorRateCutadapt = 0;
    double thresholdSW = 0;
    int lengthWindowsSW = 0;
    int seedMismatchesTrimmomatic = 0;
    int palindromeClipThresholdTrimmomatic = 0;
    int simpleClipThreshold = 0;
    int addIndexOutlier = 0;

    File samFile = new File("");
    File fastqFile = new File("");
    File fastqOutputFile = new File("");
    File adaptorFile = new File("");
    File workDir = new File("");

    try {

      if (arguments.contains("-help") || arguments.contains("-h")) {
        help(options);
      }

      // parse the command line arguments
      final CommandLine line = parser.parse(options,
          arguments.toArray(new String[arguments.size()]), true);

      // Display help
      if (line.getArgs().length == 0) {
        System.out.println(
            "ERROR : Enter the five obligatory arguments of the trim module !\n\n");
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

      if (line.hasOption("addIndexOutlier")) {
        addIndexOutlier = Integer
            .parseInt(line.getOptionValue("addIndexOutlier").toLowerCase());
      }

      if (line.hasOption("errorRateCutadapt")) {
        errorRateCutadapt = Long
            .parseLong(line.getOptionValue("errorRateCutadapt").toLowerCase());
      }

      if (line.hasOption("thresholdSW")) {
        thresholdSW =
            Long.parseLong(line.getOptionValue("thresholdSW").toLowerCase());
      }

      if (line.hasOption("lengthWindowsSW")) {
        lengthWindowsSW = Integer
            .parseInt(line.getOptionValue("lengthWindowsSW").toLowerCase());
      }

      if (line.hasOption("seedMismatchesTrimmomatic")) {
        seedMismatchesTrimmomatic = Integer.parseInt(
            line.getOptionValue("seedMismatchesTrimmomatic").toLowerCase());
      }

      if (line.hasOption("palindromeClipThresholdTrimmomatic")) {
        palindromeClipThresholdTrimmomatic = Integer
            .parseInt(line.getOptionValue("palindromeClipThresholdTrimmomatic")
                .toLowerCase());
      }

      if (line.hasOption("simpleClipThreshold")) {
        simpleClipThreshold = Integer
            .parseInt(line.getOptionValue("simpleClipThreshold").toLowerCase());
      }

      {
        String[] remainder = line.getArgs();
        if (remainder.length >= 4) {
          samFile = new File(remainder[0]);
          fastqFile = new File(remainder[1]);
          fastqOutputFile = new File(remainder[2]);
          adaptorFile = new File(remainder[3]);
          workDir = new File(remainder[4]);
        } else {
          System.out.println(
              "ERROR : Enter the five obligatory arguments of the trim module !\n\n");
          help(options);
        }
      }

    } catch (ParseException e) {
      System.out.println(
          "Error while parsing command line arguments: " + e.getMessage());
    }
    // Execute program in local mode
    run(trimmer, mode, stats, addIndexOutlier, errorRateCutadapt, thresholdSW,
        lengthWindowsSW, seedMismatchesTrimmomatic,
        palindromeClipThresholdTrimmomatic, simpleClipThreshold, samFile,
        fastqFile, fastqOutputFile, adaptorFile, workDir);
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
        .withDescription("display help").create("help"));

    options.addOption(OptionBuilder.withArgName("h").hasArg()
        .withDescription("display help").create("help"));

    options.addOption(OptionBuilder.withArgName("trimmer").hasArg()
        .withDescription(
            "name of trimmer use [cutadapt | trimmomatic] (default : cutadapt)")
        .create("trimmer"));

    options.addOption(OptionBuilder.withArgName("mode").hasArg()
        .withDescription("mode of trimming use [P | SW] (default : P)")
        .create("mode"));
    options.addOption(OptionBuilder.withArgName("stats").hasArg()
        .withDescription(
            "make somes stats on the trimming [true | false] (default : false)")
        .create("stats"));

    options.addOption(OptionBuilder.withArgName("addIndexOutlier").hasArg()
        .withDescription(
            "add more bases in addition to the outlier for P mode (default: 15")
        .create("addIndexOutlier"));

    options.addOption(OptionBuilder.withArgName("errorRateCutadapt").hasArg()
        .withDescription("error rate for cutadapt (default: 0.5")
        .create("errorRateCutadapt"));
    options.addOption(OptionBuilder.withArgName("thresholdSW").hasArg()
        .withDescription("threshold for Side-Windows processus (default: 0.8)")
        .create("thresholdSW"));
    options.addOption(OptionBuilder.withArgName("seedMismatchesTrimmomatic")
        .hasArg()
        .withDescription("seed mismatches option for Trimmomatic (default: 17)")
        .create("seedMismatchesTrimmomatic"));
    options.addOption(
        OptionBuilder.withArgName("palindromeClipThresholdTrimmomatic").hasArg()
            .withDescription(
                "palindrome clip threshold option for Trimmomatic (default: 30)")
            .create("palindromeClipThresholdTrimmomatic"));
    options.addOption(OptionBuilder.withArgName("simpleClipThreshold").hasArg()
        .withDescription(
            "simple clip threshold option for Trimmomatic (default : 7)")
        .create("simpleClipThreshold"));

    options.addOption(OptionBuilder.withArgName("samFile").hasArg()
        .withDescription("the path to the .sam file").create("samFile"));

    options.addOption(OptionBuilder.withArgName("fastqFile").hasArg()
        .withDescription("the path to the .fastq file").create("fastqFile"));

    options.addOption(OptionBuilder.withArgName("fastqOutputFile").hasArg()
        .withDescription("the path to the .fastq file")
        .create("fastqOutputFile"));

    options.addOption(OptionBuilder.withArgName("adaptorFile").hasArg()
        .withDescription("the path to the adaptor file").create("adaptorFile"));

    options.addOption(OptionBuilder.withArgName("workDir").hasArg()
        .withDescription("the work directory").create("workDir"));

    return options;
  }

  /**
   * Show command line help.
   * @param options Options of the software
   */
  private static void help(final Options options) {

    // Show help message
    final HelpFormatter formatter = new HelpFormatter();
    formatter.printHelp(Globals.APP_NAME_LOWER_CASE
        + ".sh " + ACTION_NAME + " [options] arguments \n", options);

    System.exit(0);
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
  private static void run(final String trimmer, final String mode,
      final String stats, final int addIndexOutlier,
      final double errorRateCutadapt, final double thresholdSW,
      final int lengthWindowsSW, final int seedMismatchesTrimmomatic,
      final int palindromeClipThresholdTrimmomatic,
      final int simpleClipThreshold, final File samFile, final File fastqFile,
      final File fastqOutputFile, final File adaptorFile, final File workDir) {

    try {

      getLogger().info("Sam File : " + samFile);
      getLogger().info("Fastq File: " + fastqFile);
      getLogger().info("Fastq Trimmed Output File: " + fastqOutputFile);
      getLogger().info("Adaptor File: " + adaptorFile);
      getLogger().info("Work Directory: " + workDir);

      TrimFastq trim = new TrimFastq(samFile, fastqFile, adaptorFile,
          fastqOutputFile, workDir);

      if (trimmer.contains("trimmomatic")) {
        trim.setProcessTrimmomatic(true);
      }

      if (mode.contains("SW")) {
        trim.setProcessSWTrim(true);
      }
      if (stats.contains("true")) {
        trim.setProcessStats(true);
      }
      if (thresholdSW != 0) {
        trim.setThresholdSW(thresholdSW);
      }

      if (lengthWindowsSW != 0) {
        trim.setLengthWindowSW(lengthWindowsSW);
      }

      if (addIndexOutlier != 0) {
        trim.setAddIndexOutlier(addIndexOutlier);
      }

      if (errorRateCutadapt != 0) {
        trim.setErrorRateCutadapt(errorRateCutadapt);
      }

      if (seedMismatchesTrimmomatic != 0) {
        trim.setSeedMismatchesTrimmomatic(seedMismatchesTrimmomatic);
      }

      if (palindromeClipThresholdTrimmomatic != 0) {
        trim.setPalindromeClipThresholdTrimmomatic(
            palindromeClipThresholdTrimmomatic);
      }

      if (simpleClipThreshold != 0) {
        trim.setSimpleClipThreshold(simpleClipThreshold);
      }

      trim.execution();

    } catch (Exception e3) {
      e3.printStackTrace();
    }
  }
}