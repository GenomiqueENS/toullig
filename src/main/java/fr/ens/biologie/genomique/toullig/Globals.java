package fr.ens.biologie.genomique.toullig;

import java.io.IOException;
import java.net.URL;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.Set;
import java.util.jar.Attributes;
import java.util.jar.Manifest;

public class Globals {

  private static Attributes manifestAttributes;
  private static final String MANIFEST_FILE = "/META-INF/MANIFEST.MF";

  /** The name of the application. */
  public static final String APP_NAME = "toullig";

  /** The version of the application. */
  public static final String APP_VERSION_STRING = getVersion();

  /** The build date of the application. */
  public static final String APP_BUILD_DATE = getBuiltDate();

  /** The built number of the application. */
  public static final String APP_BUILD_NUMBER = getBuiltNumber();

  /** The built commit of the application. */
  public static final String APP_BUILD_COMMIT = getBuiltCommit();

  /** The build year of the application. */
  public static final String APP_BUILD_YEAR = getBuiltYear();

  /** The name of the application. */
  public static final String APP_NAME_LOWER_CASE = APP_NAME.toLowerCase();

  /** The built host of the application. */
  public static final String APP_BUILD_HOST = getBuiltHost();

  /** The date of the copyright. */
  private static final String COPYRIGHT_DATE = "2010-" + APP_BUILD_YEAR;

  /** Platforms where the application is available. */
  public static final Set<String> AVAILABLE_BINARY_ARCH =
      Collections.unmodifiableSet(new HashSet<>(
          Arrays.asList(new String[] {"linux\tamd64", "linux\tx86_64"})));

  /** Help text. */
  public static final String HELP_TXT = Globals.APP_NAME_LOWER_CASE
      + ".sh tool [options_tool] arguments_tool \n\n"
      + "Toullig have 2 tools : \n"
      + "\t\t - fast5tofastq : Tool for read Fast5 files of minION and create the fastq.\n"
      + "\t\t - trim : Tool for trim adaptor in the fasqt of ONT.\n\n";

  /** Licence text. */
  public static final String LICENSE_TXT =
      "This program is developed under the GNU Lesser General Public License"
          + " version 2.1 or later and CeCILL-C.";

  /** About string, plain text version. */
  public static final String ABOUT_TXT = Globals.APP_NAME
      + " version " + Globals.APP_VERSION_STRING + " (" + APP_BUILD_COMMIT
      + ", " + Globals.APP_BUILD_NUMBER + ")"
      + " is a pipeline for NGS analysis.\n" + "This version has been built on "
      + APP_BUILD_DATE + ".\n\n" + "Authors:\n"
      + "  Laurent Jourdren <jourdren@biologie.ens.fr>\n"
      + "  Aurélien Birer <birer@biologie.ens.fr>\n" + "\n" + "Copyright "
      + COPYRIGHT_DATE + " IBENS genomic platform\n" + LICENSE_TXT + "\n";

  /** The welcome message. */
  public static final String WELCOME_MSG = Globals.APP_NAME
      + " version " + Globals.APP_VERSION_STRING + " (" + APP_BUILD_COMMIT
      + ", " + Globals.APP_BUILD_NUMBER + " build on " + APP_BUILD_HOST + ", "
      + Globals.APP_BUILD_DATE + ")";

  //
  // Private constants
  //

  private static final String UNKNOWN_VERSION = "UNKNOWN_VERSION";
  private static final String UNKNOWN_BUILD = "UNKNOWN_BUILD";
  private static final String UNKNOWN_DATE = "UNKNOWN_DATE";
  private static final String UNKNOWN_YEAR = "UNKNOWN_YEAR";
  private static final String UNKNOWN_BUILD_COMMIT = "UNKNOWN_COMMIT";
  private static final String UNKNOWN_BUILD_HOST = "UNKNOWN_HOST";

  //
  // Methods
  //

  private static String getVersion() {

    final String version = getManifestProperty("Specification-Version");

    return version != null ? version : UNKNOWN_VERSION;
  }

  private static String getBuiltNumber() {

    final String builtNumber = getManifestProperty("Implementation-Version");

    return builtNumber != null ? builtNumber : UNKNOWN_BUILD;
  }

  private static String getBuiltDate() {

    final String builtDate = getManifestProperty("Built-Date");

    return builtDate != null ? builtDate : UNKNOWN_DATE;
  }

  private static String getBuiltYear() {

    final String builtYear = getManifestProperty("Built-Year");

    return builtYear != null ? builtYear : UNKNOWN_YEAR;
  }

  private static String getBuiltCommit() {

    final String buildCommit = getManifestProperty("Built-Commit");

    return buildCommit != null ? buildCommit : UNKNOWN_BUILD_COMMIT;
  }

  private static String getBuiltHost() {

    final String buildHost = getManifestProperty("Built-Host");

    return buildHost != null ? buildHost : UNKNOWN_BUILD_HOST;
  }

  private static String getManifestProperty(final String propertyKey) {

    if (propertyKey == null) {
      return null;
    }

    readManifest();

    if (manifestAttributes == null) {
      return null;
    }

    return manifestAttributes.getValue(propertyKey);
  }

  private static synchronized void readManifest() {

    if (manifestAttributes != null) {
      return;
    }

    try {

      Class<?> clazz = Globals.class;
      String className = clazz.getSimpleName() + ".class";
      String classPath = clazz.getResource(className).toString();

      final String manifestPath;
      if (!classPath.startsWith("jar")) {
        // Class not from JAR

        String basePath = classPath.substring(0,
            classPath.length() - clazz.getName().length() - ".class".length());
        manifestPath = basePath + MANIFEST_FILE;

      } else {
        manifestPath = classPath.substring(0, classPath.lastIndexOf("!") + 1)
            + MANIFEST_FILE;
      }

      Manifest manifest = new Manifest(new URL(manifestPath).openStream());
      manifestAttributes = manifest.getMainAttributes();

    } catch (IOException e) {
      e.printStackTrace();
    }
  }

}
