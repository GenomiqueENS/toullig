package fr.ens.biologie.genomique.toullig;

import fr.ens.biologie.genomique.eoulsan.Common;
import fr.ens.biologie.genomique.toullig.actions.Fast5tofastqAction;
import fr.ens.biologie.genomique.toullig.actions.TrimAction;
import org.apache.commons.cli.*;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;

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
   * @return the number of options argument in the command line
   */
  public static void main(String[] args) {

    String mode = args[0];

    switch (mode) {

    case "Fast5tofastq":
      new Fast5tofastqAction().action(
          new ArrayList<String>(Arrays.asList(args)).subList(1, args.length));
      break;

    default:
      help(args);
      break;

    case "Trim":
      new TrimAction().action(
          new ArrayList<String>(Arrays.asList(args)).subList(1, args.length));
      break;

    }
  }

  /**
   * Show command line help.
   * @param args, String[] of the software
   */
  private static void help(final String[] args) {

    // Show help message
    System.out.println(
        Globals.APP_NAME_LOWER_CASE + ".sh tool [options tool] " + args);
    System.out.println("Toullig have 2 tools : \n"
        + "\t\t - Fast5tofastq : Tool for read Fast5 files of minION and create the fastq.\n"
        + "\t\t - Trim : Tool for trim adaptor in the fasqt of ONT.\n");

    Common.exit(0);
  }
}