package fr.ens.biologie.genomique.nanopore;

import java.io.File;

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
            "set a status of the fast5 file [pass|fail|failbarcode|passbarcode] to process;(default : pass)")
        .create("status"));

    options.addOption(OptionBuilder.withArgName("type").hasArg()
        .withDescription(
            "set a type of sequence [template|complement|hairpin|barcode] to process;(default : template,complement)")
        .create("type"));

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

    String status = "pass";
    String type = "template,complement";
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
        // help(options);
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
      {
        String[] remainder = line.getArgs();
        dirFast5 = new File(remainder[0]);
        dirOutputFastq = new File(remainder[1]);
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

    // FAD22491_20161011 sans barcode R9
    // FAD22487_20160810 avec barcode R9
    // FAA105486_20160617 sans barcode R7
    // File rootFast5Dir = new File("/mnt/hardMinion/FAA105486_20160617");
    //
    // File fastqDir =
    // new File("/home/birer/Bureau/nanoporetools/src/test/java/files/");

    try {
      Fast5toFastq if5 = new Fast5toFastq(dirFast5, dirOutputFastq);

      // If you want to specify the group of the status of the read in the
      // output file(ex : .._fail_complement.fastq)
      if5.setMergeAllStatusFast5(merge);
      // If the Experimental protocol is not barcoded
      if (status.contains("fail")) {
        if5.setProcessFail(true);
      }
      if (status.contains("pass")) {
        if5.setProcessPass(true);
      }
      // If the Experimental protocol is barcoded
      if (status.contains("failbarcode")) {
        if5.setProcessFailBarcode(true);
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
      if (type.contains("hairpin")) {
        if5.setSaveHairpinSequence(true);
      }
      if (type.contains("barcode")) {
        if5.setSaveBarcodeSequence(true);
      }

      if5.execution();

      try {
        LogFast5toFastq logIf5 = new LogFast5toFastq(if5, dirOutputFastq);
        logIf5.createLogConversionFastq();
      } catch (Exception e1) {
        e1.printStackTrace();
      }

    } catch (Exception e2) {
      e2.printStackTrace();
    }
  }
}