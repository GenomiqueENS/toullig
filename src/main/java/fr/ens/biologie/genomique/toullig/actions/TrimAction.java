package fr.ens.biologie.genomique.toullig.actions;

import static fr.ens.biologie.genomique.eoulsan.EoulsanLogger.getLogger;

import java.io.File;
import java.util.List;

import org.apache.commons.cli.*;

import fr.ens.biologie.genomique.eoulsan.actions.AbstractAction;
import fr.ens.biologie.genomique.toullig.Globals;
import fr.ens.biologie.genomique.toullig.trimming.TrimFastq;

/**
 * This class define the Local exec Action.
 * @since 1.0
 * @author Aurélien Birer
 */
public class TrimAction extends AbstractAction {

  /** Name of this action. */
  private static final String ACTION_NAME = "Trim";

  //
  // Action methods
  //

  /**
   * Method of TrimAction class to get the name of the action.
   * @return , a string
   */
  @Override
  public String getName() {
    return ACTION_NAME;
  }

  /**
   * Method of TrimAction class to get the desription of the action.
   * @return , a string
   */
  @Override
  public String getDescription() {
    return "execute Trim module of" + Globals.APP_NAME + " in local mode.";
  }

  /**
   * Method of TrimAction class to make the action.
   */
  @Override
  public void action(final List<String> arguments) {

    // options of the command line
    final Options options = makeOptions();

    // parser of the command line
    final CommandLineParser parser = new GnuParser();

    String trimmer = "";
    String mode = "";
    String stats = "";
    double errorRateCutadapt = 0.0;
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
      // Display help
      if (arguments.contains("-help") || arguments.contains("-h")) {
        help(options);
      }

      // Parse the command line arguments
      final CommandLine line = parser.parse(options,
          arguments.toArray(new String[arguments.size()]), true);

      // Display help if no arguments
      if (line.getArgs().length == 0) {
        System.out.println(
            "ERROR:  No argument! Please enter the five obligatory arguments of the trim module!\n\n");
        help(options);
      }

      // Get trimmer
      if (line.hasOption("trimmer")) {
        trimmer = line.getOptionValue("trimmer").toLowerCase();
      }

      // Get mode
      if (line.hasOption("mode")) {
        mode = line.getOptionValue("mode").toLowerCase();
      }

      // Get stats
      if (line.hasOption("stats")) {
        stats = line.getOptionValue("stats").toLowerCase();
      }

      // Get addIndexOutlier
      if (line.hasOption("addIndexOutlier")) {
        addIndexOutlier = Integer
            .parseInt(line.getOptionValue("addIndexOutlier").toLowerCase());
      }

      // Get errorRateCutadapt
      if (line.hasOption("errorRateCutadapt")) {
        errorRateCutadapt = Double.parseDouble(
            line.getOptionValue("errorRateCutadapt").toLowerCase());
      }

      // Get thresholdSW
      if (line.hasOption("thresholdSW")) {
        thresholdSW =
            Long.parseLong(line.getOptionValue("thresholdSW").toLowerCase());
      }

      // Get lengthWindowsSW
      if (line.hasOption("lengthWindowsSW")) {
        lengthWindowsSW = Integer
            .parseInt(line.getOptionValue("lengthWindowsSW").toLowerCase());
      }

      // Get seedMismatchesTrimmomatic
      if (line.hasOption("seedMismatchesTrimmomatic")) {
        seedMismatchesTrimmomatic = Integer.parseInt(
            line.getOptionValue("seedMismatchesTrimmomatic").toLowerCase());
      }

      // Get palindromeClipThresholdTrimmomatic
      if (line.hasOption("palindromeClipThresholdTrimmomatic")) {
        palindromeClipThresholdTrimmomatic = Integer
            .parseInt(line.getOptionValue("palindromeClipThresholdTrimmomatic")
                .toLowerCase());
      }

      // Get simpleClipThreshold
      if (line.hasOption("simpleClipThreshold")) {
        simpleClipThreshold = Integer
            .parseInt(line.getOptionValue("simpleClipThreshold").toLowerCase());
      }

      // Get arguments
      {
        String[] remainder = line.getArgs();
        if (remainder.length >= 4) {

          // Get sam File
          samFile = new File(remainder[0]);

          // Get fastq File
          fastqFile = new File(remainder[1]);

          // Get fastq Output File
          fastqOutputFile = new File(remainder[2]);

          // Get adaptor File
          adaptorFile = new File(remainder[3]);

          // Get working Directory
          workDir = new File(remainder[4]);

        } else {
          System.out.println(
              "ERROR: Enter the five obligatory arguments of the trim module!\n\n");

          // display help
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

    // add option for help
    options.addOption(OptionBuilder.withArgName("help").hasArg()
        .withDescription("display help").create("help"));

    // add option for help
    options.addOption(OptionBuilder.withArgName("h").hasArg()
        .withDescription("display help").create("help"));

    // add option for trimmer
    options.addOption(OptionBuilder.withArgName("trimmer").hasArg()
        .withDescription(
            "name of trimmer use [cutadapt | trimmomatic | no] (default : cutadapt)")
        .create("trimmer"));

    // add option for mode
    options.addOption(OptionBuilder.withArgName("mode").hasArg()
        .withDescription("mode of cutadaptTrimming use [P | SW] (default : P)")
        .create("mode"));

    // add option for stats
    options.addOption(OptionBuilder.withArgName("stats").hasArg()
        .withDescription(
            "make somes stats on the cutadaptTrimming [true | false] (default : false)")
        .create("stats"));

    // add option for add Index to the Outlier during the trim
    options.addOption(OptionBuilder.withArgName("addIndexOutlier").hasArg()
        .withDescription(
            "add more bases in addition to the outlier for P mode (default: 15")
        .create("addIndexOutlier"));

    // add option for error Rate Cutadapt
    options.addOption(OptionBuilder.withArgName("errorRateCutadapt").hasArg()
        .withDescription("error rate for cutadapt (default: 0.5")
        .create("errorRateCutadapt"));

    // add option for threshold Side Window
    options.addOption(OptionBuilder.withArgName("thresholdSW").hasArg()
        .withDescription("threshold for Side-Windows processus (default: 0.8)")
        .create("thresholdSW"));

    // add option for seed mismatches for trimmomatic
    options.addOption(OptionBuilder.withArgName("seedMismatchesTrimmomatic")
        .hasArg()
        .withDescription("seed mismatches option for Trimmomatic (default: 17)")
        .create("seedMismatchesTrimmomatic"));

    // add option for palindrome Clip Threshold for trimmomatic
    options.addOption(
        OptionBuilder.withArgName("palindromeClipThresholdTrimmomatic").hasArg()
            .withDescription(
                "palindrome clip threshold option for Trimmomatic (default: 30)")
            .create("palindromeClipThresholdTrimmomatic"));

    // add option for simple Clip Threshold for trimmomatic
    options.addOption(OptionBuilder.withArgName("simpleClipThreshold").hasArg()
        .withDescription(
            "simple clip threshold option for Trimmomatic (default : 7)")
        .create("simpleClipThreshold"));

    // return options
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
            + ".sh " + ACTION_NAME
            + " [options] SAM_FILE FASTQ_FILE FASTQ_OUTPUT_FILE ADAPTOR_FILE WORK_DIR \n",
        options);

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

      // Logger of the action
      getLogger().info("Sam File : " + samFile);
      getLogger().info("Fastq File: " + fastqFile);
      getLogger().info("Fastq Trimmed Output File: " + fastqOutputFile);
      getLogger().info("Adaptor File: " + adaptorFile);
      getLogger().info("Work Directory: " + workDir);

      // Call the constructor with the arguments
      TrimFastq trim = new TrimFastq(samFile, fastqFile, adaptorFile,
          fastqOutputFile, workDir, trimmer, mode);

      // if the stats on the trim will be display
      if (stats.contains("true")) {
        trim.setProcessStatsCutadapt();
      }

      // set the threshold for Side Window method
      if (thresholdSW != 0) {
        trim.setThresholdSideWindow(thresholdSW);
      }

      // set the length of the window for Side Window method
      if (lengthWindowsSW != 0) {
        trim.setLengthWindowSideWindow(lengthWindowsSW);
      }

      // set the number of add index to the outlier
      trim.setAddIndexOutlier(addIndexOutlier);

      // set the error rate fort cutadapt trimmmer
      if (errorRateCutadapt != 0) {
        trim.setErrorRateCutadapt(errorRateCutadapt);
      }

      // set the seed mismatches for trimmomatic trimmer
      if (seedMismatchesTrimmomatic != 0) {
        trim.setSeedMismatchesTrimmomatic(seedMismatchesTrimmomatic);
      }

      // set the palindrome clip Threshold for trimmomatic trimmer
      if (palindromeClipThresholdTrimmomatic != 0) {
        trim.setPalindromeClipThresholdTrimmomatic(
            palindromeClipThresholdTrimmomatic);
      }

      // set the simple Clip Threshold for trimmomatic trimmer
      if (simpleClipThreshold != 0) {
        trim.setSimpleClipThreshold(simpleClipThreshold);
      }

      // execute the trim
      trim.execution();

    } catch (Exception e3) {
      e3.printStackTrace();
    }
  }
}
