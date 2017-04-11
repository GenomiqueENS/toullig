package fr.ens.biologie.genomique.toullig.actions;

import static fr.ens.biologie.genomique.eoulsan.EoulsanLogger.getLogger;

import java.io.File;
import java.util.Date;
import java.util.List;

import org.apache.commons.cli.*;

import fr.ens.biologie.genomique.eoulsan.actions.AbstractAction;
import fr.ens.biologie.genomique.toullig.fast5tofastq.Fast5ToFastq;
import fr.ens.biologie.genomique.toullig.fast5tofastq.Fast5ToFastqLogger;
import fr.ens.biologie.genomique.toullig.Globals;

/**
 * This class define the Local exec Action.
 * @since 1.0
 * @author Aurélien Birer
 */
public class Fast5tofastqAction extends AbstractAction {

  /** Name of this action. */
  private static final String ACTION_NAME = "Fast5tofastq";

  //
  // Action methods
  //

  /**
   * Method of Fast5tofastqAction class to get the name of the action.
   * @return , a string
   */
  @Override
  public String getName() {
    return ACTION_NAME;
  }

  /**
   * Method of Fast5tofastqAction class to get the desription of the action.
   * @return , a string
   */
  @Override
  public String getDescription() {
    return "execute Fast5tofastq module of "
        + Globals.APP_NAME + " in local mode.";
  }

  /**
   * Method of Fast5tofastqAction class to make the action.
   */
  @Override
  public void action(final List<String> arguments) {

    // options of the command line
    final Options options = makeOptions();

    // parser of the command line
    final CommandLineParser parser = new GnuParser();

    String status = "pass";
    String type = "transcript";
    String compress = "";
    File dirFast5 = null;
    File dirOutputFastq = null;
    boolean merge = false;

    try {

      // Display help
      if (arguments.contains("-help") || arguments.contains("-h")) {
        help(options);
      }

      // parse the command line arguments
      final CommandLine line = parser.parse(options,
          arguments.toArray(new String[arguments.size()]), true);

      // Display help if no arguments
      if (line.getArgs().length == 0) {
        System.out.println(
            "ERROR: No argument! Please enter the two obligatory arguments\n\n");
        help(options);
      }

      // Get status
      if (line.hasOption("status")) {
        status = line.getOptionValue("status").toLowerCase();
      }

      // Get type
      if (line.hasOption("type")) {
        type = line.getOptionValue("type").toLowerCase();
      }
      // Get compression format
      if (line.hasOption("compress")) {
        compress = line.getOptionValue("compress").toLowerCase();
      }
      // Get arguments
      {
        String[] remainder = line.getArgs();
        if (remainder.length >= 2) {

          // Get directory of Fast5 run
          dirFast5 = new File(remainder[0]);

          // Get directory of Output Fastq
          dirOutputFastq = new File(remainder[1]);

        } else {
          System.out.println(
              "ERROR: Enter the two obligatory arguments of the directory of the basecalled run and the output directory for fastq!\n\n");
          // display help
          help(options);
        }
      }

      // Get mergeSequence options
      if (line.hasOption("mergeSequence")) {
        merge = Boolean.parseBoolean(line.getOptionValue("merge"));

      }

    } catch (ParseException e) {
      System.out.println(
          "Error while parsing command line arguments: " + e.getMessage());
    }

    // Execute program in local mode
    run(status, type, compress, dirFast5, dirOutputFastq, merge, arguments);
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

    // add option for status
    options.addOption(OptionBuilder.withArgName("status").hasArg()
        .withDescription(
            "set a status of the fast5 file [pass|fail|unclassified] to process;(default: none)")
        .create("status"));

    // add option for type
    options.addOption(OptionBuilder.withArgName("type").hasArg()
        .withDescription(
            "set a type of sequence [template|complement|consensus|transcript] to obtain;(default: transcript)")
        .create("type"));

    // add option for compress
    options.addOption(OptionBuilder.withArgName("compress").hasArg()
        .withDescription(
            "set a compression for the output fastq [gzip|bzip2] (default: none)")
        .create("compress"));

    // add option for mergeSequence
    options.addOption(OptionBuilder.withArgName("mergeSequence").hasArg()
        .withDescription(
            "merge the sequence of status choose [true/false];(default: false)")
        .create("merge"));

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
            + ".sh " + ACTION_NAME + "[options] FAST5_DIR FASTQ_OUTPUT_DIR\n",
        options);

    System.exit(0);
  }

  //
  // Execution
  //

  /**
   * Run Toullig Fast5tofastq
   * @param status, the status of the fast5
   * @param type, the type of the sequence
   * @param compress, the type of compression
   * @param dirFast5, the root directory of the run fast5
   * @param dirOutputFastq, the output directory for fastq
   * @param merge, boolean for assemble all data
   */
  private static void run(final String status, final String type,
      final String compress, final File dirFast5, final File dirOutputFastq,
      final boolean merge, List<String> arguments) {

    // Get the Begin Date of the action
    Date beginDate = new Date();
    try {

      // Call the constructor with the arguments
      Fast5ToFastq if5 = new Fast5ToFastq(dirFast5, dirOutputFastq);

      // If you want to specify the group of the status of the read in the
      // output file(ex : .._fail_complement.fastq)
      if5.setMergeAllStatusFast5(merge);

      // If the Experimental protocol is not barcoded
      if (status.contains("fail")) {
        if5.setProcessFail();
      }

      // set the pass fast5 to process
      if (status.contains("pass")) {
        if5.setProcessPass();
      }
      // If the Experimental protocol is barcoded
      if (status.contains("unclassified")) {
        if5.setProcessUnclassified();
      }

      // set Compress format bzip2
      if (type.contains("template")) {
        if5.setSaveTemplateSequence();
      }

      // set Compress format bzip2
      if (type.contains("complement")) {
        if5.setSaveComplementSequence();
      }

      // set Compress format bzip2
      if (type.contains("consensus")) {
        if5.setSaveConsensusSequence();
      }

      // set Compress format bzip2
      if (type.contains("transcript")) {
        if5.setSaveTranscriptSequence();
      }
      // set Compress format gzip
      if (compress.contains("gzip")) {
        if5.setGzipCompression();
      }

      // set Compress format bzip2
      if (compress.contains("bzip2")) {
        if5.setBZip2Compression();
      }

      // Logger of the action
      getLogger().info("Fast5 Run Directory: " + dirFast5);
      getLogger().info("Fastq Output Directory: " + dirOutputFastq);

      // Execution to the read of fast5 to fastq
      if5.execute();

      // Get the end Date of the execution
      Date endDate = new Date();

      // Write of few logs files
      try {

        Fast5ToFastqLogger logIf5 = new Fast5ToFastqLogger(if5, dirOutputFastq);
        logIf5.createLogConversionFastq(beginDate, endDate, arguments);
        logIf5.createLogCorruptFile();
        logIf5.createLogWorkflow();

      } catch (Exception e1) {
        e1.printStackTrace();
      }

    } catch (Exception e2) {
      e2.printStackTrace();
    }
  }
}