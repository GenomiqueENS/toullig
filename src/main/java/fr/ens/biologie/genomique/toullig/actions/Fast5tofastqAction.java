package fr.ens.biologie.genomique.toullig.actions;

import static fr.ens.biologie.genomique.eoulsan.EoulsanLogger.getLogger;

import java.io.File;
import java.util.Date;
import java.util.List;

import fr.ens.biologie.genomique.eoulsan.actions.AbstractAction;
import fr.ens.biologie.genomique.toullig.Fast5ToFastq;
import fr.ens.biologie.genomique.toullig.Fast5ToFastqLogger;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import fr.ens.biologie.genomique.eoulsan.Common;
import fr.ens.biologie.genomique.toullig.Globals;

/**
 * This class define the Local exec Action.
 * @since 1.0
 * @author Aurélien Birer
 */
public class Fast5tofastqAction extends AbstractAction {

  /** Name of this action. */
  public static final String ACTION_NAME = "Fast5tofastq";

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

    String status = "pass";
    String type = "transcript";
    String compress = "";
    File dirFast5 = null;
    File dirOutputFastq = null;
    boolean merge = false;

    try {

      // parse the command line arguments
      final CommandLine line = parser.parse(options,
          arguments.toArray(new String[arguments.size()]), true);

      // Display help
      if (line.hasOption("help")) {
        help(options);
      }

      // Set status
      if (line.hasOption("status")) {
        status = line.getOptionValue("status").toLowerCase();
      }
      // Set type
      if (line.hasOption("type")) {
        type = line.getOptionValue("type").toLowerCase();
      }
      // Set compression format
      if (line.hasOption("compress")) {
        compress = line.getOptionValue("compress").toLowerCase();
      }
      {
        String[] remainder = line.getArgs();
        if (remainder.length >= 2) {
          dirFast5 = new File(remainder[0]);
          dirOutputFastq = new File(remainder[1]);
        }
      }

      // Set the root directory of the Fast5 run
      if (line.hasOption("mergeSequence")) {
        merge = Boolean.parseBoolean(line.getOptionValue("merge"));

      }

    } catch (ParseException e) {
      Common.errorExit(e,
          "Error while parsing command line arguments: " + e.getMessage());
    }

    // Execute program in local mode
    run(status, type, compress, dirFast5, dirOutputFastq, merge);
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

    options.addOption(OptionBuilder.withArgName("status").hasArg()
        .withDescription(
            "set a status of the fast5 file [pass|fail|unclassified|passbarcode] to process;(default : pass)")
        .create("status"));

    options.addOption(OptionBuilder.withArgName("type").hasArg()
        .withDescription(
            "set a type of sequence [template|complement|consensus|transcript] to obtain ;(default : transcript)")
        .create("type"));

    options.addOption(OptionBuilder.withArgName("compress").hasArg()
        .withDescription("set a compression for the output fastq [GZIP|BZIP2]")
        .create("compress"));

    options
        .addOption(OptionBuilder.withArgName("rootDirectoryFast5run").hasArg()
            .withDescription("the path to the main directory of the minION run")
            .create("dirFast5"));

    options.addOption(OptionBuilder.withArgName("outputDirectoryFastq").hasArg()
        .withDescription("the path to the output of Fastq sequences")
        .create("dirOutputFastq"));

    options.addOption(OptionBuilder.withArgName("mergeSequence").hasArg()
        .withDescription(
            "merge the sequence of status choose [true/false];(default : false)")
        .create("merge"));

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
      final boolean merge) {

    Date dateDeb = new Date();
    try {
      Fast5ToFastq if5 = new Fast5ToFastq(dirFast5, dirOutputFastq);

      // If you want to specify the group of the status of the read in the
      // output file(ex : .._fail_complement.fastq)
      if5.setMergeAllStatusFast5(merge);
      // If the Experimental protocol is not barcoded
      if (status.contains("fail")) {
        if5.setProcessFail(true);
      }
      if (status.contains("pass") && !status.contains("passbarcode")) {
        if5.setProcessPass(true);
      }
      // If the Experimental protocol is barcoded
      if (status.contains("unclassified")) {
        if5.setProcessUnclassified(true);
      }
      if (status.contains("passbarcode")) {
        if5.setProcessPassBarcode(true);
      }
      // List of sequences :
      if (type.contains("template")) {
        if5.setSaveTemplateSequence(true);
      }
      if (type.contains("complement")) {
        if5.setSaveComplementSequence(true);
      }
      if (type.contains("consensus")) {
        if5.setSaveConsensusSequence(true);
      }
      if (type.contains("transcript")) {
        if5.setSaveTranscriptSequence(true);
      }
      // Compress format
      if (compress.contains("gzip")) {
        if5.setGzipCompression(true);
      }
      if (compress.contains("bzip2")) {
        if5.setBZip2Compression(true);
      }

      getLogger().info("Fast5 Run Directory : " + dirFast5);
      getLogger().info("Fastq Output Directory: " + dirOutputFastq);

      // Execution to the read of fast5 to fastq
      if5.execute();

      Date dateEnd = new Date();
      // Write of few logs files
      try {
        Fast5ToFastqLogger logIf5 = new Fast5ToFastqLogger(if5, dirOutputFastq);
        logIf5.createLogConversionFastq(dateDeb, dateEnd);
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