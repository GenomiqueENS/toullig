package fr.ens.biologie.genomique.toullig;

import java.util.ArrayList;
import java.util.Arrays;

import fr.ens.biologie.genomique.eoulsan.Common;
import fr.ens.biologie.genomique.toullig.actions.GtftogpdAction;
import org.apache.commons.cli.*;

import fr.ens.biologie.genomique.toullig.actions.Fast5tofastqAction;
import fr.ens.biologie.genomique.toullig.actions.TrimAction;

/**
 * Main class of nanoporetools
 * @author Aurelien Birer
 */
public class Main {

  /**
   * Show a message and then exit.
   * @param message the message to show
   */
  public static void showErrorMessageAndExit(final String message) {

    System.err.println(message);
    System.exit(1);
  }

  /**
   * Show a message and then exit.
   * @param message the message to show
   */
  public static void showMessageAndExit(final String message) {

    System.out.println(message);
    System.exit(0);
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

    // add option for display version
    options.addOption("version", false, "show version of the software");

    // add option for display about
    options.addOption("about", false,
        "display information about this software");

    // add option for display help
    options.addOption("h", "help", false, "display this help");

    // add option for display license
    options.addOption("license", false,
        "display information about the license of this software");

    return options;
  }

  /**
   * The main function parse the options of the command line
   */
  public static void main(String[] args) {

    // make options
    final Options options = makeOptions();

    // make parser of the command line
    final CommandLineParser parser = new GnuParser();

    try {

      // parse the command line arguments
      final CommandLine line = parser.parse(options, args, true);

      // test if no arguments is given
      if (args.length == 0) {

        showErrorMessageAndExit(
            "ERROR: Toullig need in the first argument the tool what you want use!\nSee the help with "
                + Globals.APP_NAME_LOWER_CASE + ".sh -h, --help\n");

      }

      // display help
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

      // process fast5tofastq module
      case "fast5tofastq":

        new Fast5tofastqAction().action(
            new ArrayList<>(Arrays.asList(args)).subList(1, args.length));

        break;

      // display help
      default:

        showMessageAndExit(
            "ERROR: The name of the tool is not correct!\nSee the help with "
                + Globals.APP_NAME_LOWER_CASE + ".sh -h, --help\n");

        break;

      // process trim module
      case "trim":

        new TrimAction().action(
            new ArrayList<>(Arrays.asList(args)).subList(1, args.length));

        break;

      // process gtftogpd module
        case "gtftogpd":

        new GtftogpdAction().action(
                new ArrayList<>(Arrays.asList(args)).subList(1, args.length));

        break;


      // // process test module
      // case "test":
      // new SamDiffAnnot();
      // break;

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
    System.out.println(Globals.HELP_TXT);

    System.exit(0);
  }
}