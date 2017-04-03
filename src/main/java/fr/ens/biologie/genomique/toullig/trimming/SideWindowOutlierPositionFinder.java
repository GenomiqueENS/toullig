package fr.ens.biologie.genomique.toullig.trimming;

import java.io.File;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Class to execute the SideWindow method to find outliers. Created by birer on
 * 29/03/17.
 */
public class SideWindowOutlierPositionFinder implements OutlierPositionFinder {

  private final int lengthWindowsSideWindow;
  private final double thresholdSideWindow;
  private final Map<String, String[]> workTrimmingMap;

  public SideWindowOutlierPositionFinder(int lengthWindowsSideWindow,
      double thresholdSideWindow, Map<String, String[]> workTrimmingMap) {

    this.lengthWindowsSideWindow = lengthWindowsSideWindow;
    this.thresholdSideWindow = thresholdSideWindow;
    this.workTrimmingMap = workTrimmingMap;
  }

  /**
   * Method of the class SideWindowOutlierPositionFinder to trim sequence to
   * create sequences files for cutadapt.
   * @param fastaLeftOutlierFile , a fasta File for left outlier
   * @param fastaRightOutlierFile , a fasta File for right outlier
   * @param trimmer, the trimmer interface
   */
  public void findOutliers(File fastaLeftOutlierFile,
      File fastaRightOutlierFile, Trimmer trimmer) {

    Pattern oneCodeCigarPattern = Pattern.compile("(([0-9]*[A-Z]).*)");
    Pattern lengthOneCodeCigarPattern = Pattern.compile("([0-9]*)(.)");
    int lengthSequenceCigar;

    for (String id : this.workTrimmingMap.keySet()) {

      String sequenceCigarBinary = "";
      String[] tabValue = this.workTrimmingMap.get(id);
      String sequence = tabValue[0];
      String cigar = tabValue[2];
      String quality = tabValue[1];

      if (!"*".equals(cigar)) {

        while (cigar.length() != 0) {

          Matcher oneCodeCigarMatcher = oneCodeCigarPattern.matcher(cigar);

          // if one Code Cigar is found
          if (oneCodeCigarMatcher.matches()) {

            // for each one Code Cigar
            String oneCodeCigar = oneCodeCigarMatcher.group(2);
            cigar = cigar.substring(oneCodeCigarMatcher.group(2).length(),
                cigar.length());
            Matcher lengthOneCodeCigarMatcher =
                lengthOneCodeCigarPattern.matcher(oneCodeCigar);

            // if a length of one Code Cigar is found
            if (lengthOneCodeCigarMatcher.matches()) {

              lengthSequenceCigar =
                  Integer.parseInt(lengthOneCodeCigarMatcher.group(1)) - 1;

              for (int i = 0; i <= lengthSequenceCigar; i++) {

                if (lengthOneCodeCigarMatcher.group(2).equals("M")) {
                  sequenceCigarBinary = sequenceCigarBinary + 1;
                  continue;
                }
                if (lengthOneCodeCigarMatcher.group(2).equals("N")
                    || lengthOneCodeCigarMatcher.group(2).equals("D")) {
                  continue;
                }
                if (!lengthOneCodeCigarMatcher.group(2).equals("N")
                    || !lengthOneCodeCigarMatcher.group(2).equals("M")) {
                  sequenceCigarBinary = sequenceCigarBinary + 0;
                }
              }
            }
          }
        }

        int leftLengthOutlier = sideWindowsLeft(sequenceCigarBinary);
        int rightLengthOutlier = sideWindowsRight(sequenceCigarBinary);
        tabValue[3] = "" + leftLengthOutlier;
        tabValue[4] = "" + rightLengthOutlier;

        System.out.println(leftLengthOutlier
            + "    " + rightLengthOutlier + "     " + sequence.length());

        // pre-process trimming
        trimmer.preProcessTrimming(leftLengthOutlier, rightLengthOutlier,
            sequence, id, quality);
      }
    }

  }

  //
  // get index by side window
  //

  /**
   * Method of the class SideWindowOutlierPositionFinder to the length of the
   * left outlier with a binary CIGAR sequence.
   * @param sequenceCigarBinary, a binaire CIGAR sequence
   * @return int, the length of the outlier
   */
  private int sideWindowsLeft(String sequenceCigarBinary) {

    if (sequenceCigarBinary == null) {
      return 0;
    }
    String window;
    int length = 0;
    for (int i = 0; i <= sequenceCigarBinary.length(); i++) {
      if (i == (sequenceCigarBinary.length()
          - this.lengthWindowsSideWindow - 2)) {
        break;
      }
      window =
          sequenceCigarBinary.substring(i, i + this.lengthWindowsSideWindow);
      if (sumWindowCIGAR(window) >= this.thresholdSideWindow) {
        length = i + this.lengthWindowsSideWindow;
        break;
      }
    }
    return length;
  }

  /**
   * Method of the class SideWindowOutlierPositionFinder to the length of the
   * right outlier with a binary CIGAR sequence.
   * @param sequenceCigarBinary, a binaire CIGAR sequence
   * @return int, the length of the outlier
   */
  private int sideWindowsRight(String sequenceCigarBinary) {

    if (sequenceCigarBinary == null) {
      return 0;
    }
    String window;
    int length = 0;
    for (int i = sequenceCigarBinary.length(); i >= 0; i--) {
      if (i == this.lengthWindowsSideWindow) {
        break;
      }

      window =
          sequenceCigarBinary.substring(i - this.lengthWindowsSideWindow, i);
      if (sumWindowCIGAR(window) >= this.thresholdSideWindow) {
        if (i >= sequenceCigarBinary.length() - this.lengthWindowsSideWindow) {
          length = sequenceCigarBinary.length() - i;
        } else {
          length =
              sequenceCigarBinary.length() - i - this.lengthWindowsSideWindow;
        }
        break;
      }
    }
    return length;
  }

  /**
   * Method of the class SideWindowOutlierPositionFinder to sumWindowCIGAR a
   * sequence CIGAR encode in binary.
   * @param windows, a String of CIGAR sequence encode in binary
   * @return double, the sumWindowCIGAR of the CIGAR window
   */
  private double sumWindowCIGAR(String windows) {
    char[] arrayWindows = windows.toCharArray();
    int sum = 0;
    for (int i = 0; i <= arrayWindows.length - 1; i++) {
      if (arrayWindows[i] == 1) {
        sum += 1;
      }
    }
    return (double) sum / windows.length();
  }
}
