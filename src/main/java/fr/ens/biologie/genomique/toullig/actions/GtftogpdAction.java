package fr.ens.biologie.genomique.toullig.actions;

import static fr.ens.biologie.genomique.eoulsan.EoulsanLogger.getLogger;

import java.io.File;
import java.util.List;

import fr.ens.biologie.genomique.toullig.GtftoGdp.GtftoGpd;
import org.apache.commons.cli.*;

import fr.ens.biologie.genomique.eoulsan.actions.AbstractAction;
import fr.ens.biologie.genomique.toullig.Globals;

/**
 * This class define the Local exec Action.
 * @since 1.0
 * @author Aurélien Birer
 */
public class GtftogpdAction extends AbstractAction {

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

    File gtfFile = new File("");
    File gpdOutputFile = new File("");

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
            "ERROR:  No argument! Please enter the 2 obligatory arguments of the gtftogpd module!\n\n");
        help(options);
      }

      // Get arguments
      {
        String[] remainder = line.getArgs();
        if (remainder.length >= 1) {

          // Get sam File
            gtfFile = new File(remainder[0]);

          // Get fastq File
            gpdOutputFile = new File(remainder[1]);

        } else {
          System.out.println(
              "ERROR: Enter the 2 obligatory arguments of the gtftogpd module!\n\n");

          // display help
          help(options);
        }
      }

    } catch (ParseException e) {
      System.out.println(
          "Error while parsing command line arguments: " + e.getMessage());
    }
    // Execute program in local mode
    run(gtfFile,gpdOutputFile);
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
            + " GTF_FILE GPD_NAME_OUTPUT_FILE \n",
        options);

    System.exit(0);
  }

  //
  // Execution
  //

  /**
   * Run Toullig Trim
   * @param gtfFile, a gtf file
   * @param gpdOutputFile, a gpd file
   */
  private static void run(final File gtfFile, final File gpdOutputFile) {

    try {

      // Logger of the action
      getLogger().info("GTF File : " + gtfFile);
      getLogger().info("GPD output File: " + gpdOutputFile);

      // Call the constructor with the arguments
        GtftoGpd exec = new GtftoGpd(gtfFile,gpdOutputFile);

      // execute the translator
        exec.execution();

    } catch (Exception e) {
      e.printStackTrace();
    }
  }

}
