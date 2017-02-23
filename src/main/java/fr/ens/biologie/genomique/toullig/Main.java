package fr.ens.biologie.genomique.toullig;

import java.io.File;
import java.util.Date;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;


/**
 * Main class of nanoporetools
 * @author birer
 */
public abstract class Main {

  /**
   * Exit the application.
   * @param exitCode exit code
   */
  public static void exit(final int exitCode) {

    System.exit(exitCode);
  }

  /**
   * Show a message and then exit.
   * @param message the message to show
   */
  public static void showErrorMessageAndExit(final String message) {

    System.err.println(message);
    exit(1);
  }

  /**
   * Show a message and then exit.
   * @param message the message to show
   */
  public static void showMessageAndExit(final String message) {

    System.out.println(message);
    exit(0);
  }

  /**
   * Create options for command line
   * @return an Options object
   */
  @SuppressWarnings("static-access")
  private static Options makeOptions() {

    // create Options object
    final Options options = new Options();

    options.addOption("version", false, "show version of the software");
    options.addOption("about", false,
        "display information about this software");
    options.addOption("h", "help", false, "display this help");
    options.addOption("license", false,
        "display information about the license of this software");

    options.addOption(OptionBuilder.withArgName("status").hasArg()
        .withDescription(
            "set a status of the fast5 file [pass|fail|unclassified|passbarcode] to process;(default : pass)")
        .create("status"));

    options.addOption(OptionBuilder.withArgName("type").hasArg()
        .withDescription(
            "set a type of sequence [template|complement|consensus|transcript] to obtain ;(default : transcript)")
        .create("type"));
    
    options.addOption(OptionBuilder.withArgName("compress").hasArg()
        .withDescription(
            "set a compression for the output fastq [GZIP|BZIP2]")
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
   * The main function parse the options of the command line
   * @return the number of options argument in the command line
   */
  public static void main(String[] args) {
    
    Date dateDeb = new Date();

    String status = "pass";
    String type = "transcript";
    String compress = "";
    File dirFast5 = null;
    File dirOutputFastq = null;
    boolean merge = false;

    // create Options
    final Options options = makeOptions();

    // create CommandLineParser
    CommandLineParser parser = new GnuParser();

    try {
      // parse the command line arguments
      CommandLine line = parser.parse(options, args, true);
      // Help option
      if (line.hasOption("help")) {
        showMessageAndExit(Globals.HELP_TXT);
      }
      // About option
      if (line.hasOption("about")) {
        showMessageAndExit(Globals.ABOUT_TXT);
      }
      // Version option
      if (line.hasOption("version")) {
        showMessageAndExit(Globals.WELCOME_MSG);
      }
      // Licence option
      if (line.hasOption("license")) {
        showMessageAndExit(Globals.LICENSE_TXT);
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
      e.printStackTrace();
    }
    if (args == null) {
      showErrorMessageAndExit("This program needs one argument."
          + " Use the -h option to get more information.\n");
    }

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
      //Compress format
      if(compress.contains("gzip")) {
        if5.setGzipCompression(true);
      }
      if(compress.contains("bzip2")) {
        if5.setBZip2Compression(true);
      }
      
      //Execution to the read of fast5 to fastq

      if5.execute();
      
      Date dateEnd = new Date();
      
      //Write of few logs files
      try {
        Fast5ToFastqLogger logIf5 = new Fast5ToFastqLogger(if5, dirOutputFastq);
        logIf5.createLogConversionFastq(dateDeb,dateEnd);
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