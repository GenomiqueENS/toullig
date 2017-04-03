package fr.ens.biologie.genomique.toullig;

import java.util.ArrayList;
import java.util.Arrays;

import fr.ens.biologie.genomique.eoulsan.Common;
import org.apache.commons.cli.*;

import fr.ens.biologie.genomique.toullig.actions.Fast5tofastqAction;
import fr.ens.biologie.genomique.toullig.actions.TrimAction;

/**
 * Main class of nanoporetools
 * @author birer
 */
abstract class Main {

  /**
   * Exit the application.
   * @param exitCode exit code
   */
  private static void exit(final int exitCode) {

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
  private static Options makeOptions() {

    // create Options object
    final Options options = new Options();

    //
    // General
    //

    options.addOption("version", false, "show version of the software");
    options.addOption("about", false,
        "display information about this software");
    options.addOption("h", "help", false, "display this help");
    options.addOption("license", false,
        "display information about the license of this software");

    options.addOption(OptionBuilder.withArgName("mode").hasArg()
        .withDescription("mode of toullig [fast5tofastq|trim]").create("mode"));

    return options;
  }

  /**
   * The main function parse the options of the command line
   */
  public static void main(String[] args) {

    final Options options = makeOptions();
    final CommandLineParser parser = new GnuParser();

    try {

      // parse the command line arguments
      final CommandLine line = parser.parse(options, args, true);

      if (args.length == 0) {
        System.out.println(
            "ERROR: Toullig need in the first argument the tool what you want use!");
        System.out.println("See the help with "
            + Globals.APP_NAME_LOWER_CASE + ".sh -h, --help\n");
        exit(1);

      }

      if (line.hasOption("help")) {
        help();
      }

      // About option
      if (line.hasOption("about")) {
        Common.showMessageAndExit(Globals.ABOUT_TXT);
      }

      // Version option
      if (line.hasOption("version")) {
        Common.showMessageAndExit(Globals.WELCOME_MSG);
      }

      // Licence option
      if (line.hasOption("license")) {
        Common.showMessageAndExit(Globals.LICENSE_TXT);
      }

      String mode = args[0];

      switch (mode) {

      case "fast5tofastq":
        new Fast5tofastqAction().action(
            new ArrayList<>(Arrays.asList(args)).subList(1, args.length));
        break;

      default:
        System.out.println("ERROR: The name of the tool is not correct!");
        System.out.println("See the help with "
            + Globals.APP_NAME_LOWER_CASE + ".sh -h, --help\n");
        break;

      case "trim":
        new TrimAction().action(
            new ArrayList<>(Arrays.asList(args)).subList(1, args.length));
        break;

      }
    } catch (ParseException e) {
      e.printStackTrace();
    }
  }

  /**
   * Show command line help.
   */
  private static void help() {

    // Show help message
    System.out.println(
        Globals.APP_NAME_LOWER_CASE + ".sh tool [options_tool] arguments_tool");
    System.out.println("Toullig have 2 tools : \n"
        + "\t\t - fast5tofastq : Tool for read Fast5 files of minION and create the fastq.\n"
        + "\t\t - trim : Tool for trim adaptor in the fasqt of ONT.\n\n");

    exit(0);
  }
}