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
  private final Map<String, InformationRead> workTrimmingMap;

  /**
   * Constructor of the SideWindow class.
   * @param lengthWindowsSideWindow, a int of the length of the window for the
   *          side-window outlier postion finder
   * @param thresholdSideWindow, a double of the threshold for the side-window
   *          outlier postion finder
   * @param workTrimmingMap, a Map of working information
   */
  public SideWindowOutlierPositionFinder(int lengthWindowsSideWindow,
      double thresholdSideWindow,
      Map<String, InformationRead> workTrimmingMap) {

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

    // pattern to have one Cigar Code
    Pattern oneCodeCigarPattern = Pattern.compile("(([0-9]*[A-Z]).*)");

    // pattern to have the length of one Cigar Code
    Pattern lengthOneCodeCigarPattern = Pattern.compile("([0-9]*)(.)");

    int lengthSequenceCigar;

    // get all read of the working Map
    for (String id : this.workTrimmingMap.keySet()) {

      String sequenceCigarBinary = "";

      // get information of the read
      InformationRead informationRead = this.workTrimmingMap.get(id);
      String sequence = informationRead.sequence;
      String cigar = informationRead.cigar;
      String quality = informationRead.quality;

      // test if the read is unmapped
      if (!"*".equals(cigar)) {

        // while the length of the cigar is not null
        while (cigar.length() != 0) {

          // match the one Code Cigar
          Matcher oneCodeCigarMatcher = oneCodeCigarPattern.matcher(cigar);

          // if one Code Cigar is found
          if (oneCodeCigarMatcher.matches()) {

            // for each one Code Cigar
            String oneCodeCigar = oneCodeCigarMatcher.group(2);

            // decremente the length of the Cigar code
            cigar = cigar.substring(oneCodeCigarMatcher.group(2).length(),
                cigar.length());

            // match the length of one Code Cigar
            Matcher lengthOneCodeCigarMatcher =
                lengthOneCodeCigarPattern.matcher(oneCodeCigar);

            // if a length of one Code Cigar is found
            if (lengthOneCodeCigarMatcher.matches()) {

              // get the length of the cigar sequence
              lengthSequenceCigar =
                  Integer.parseInt(lengthOneCodeCigarMatcher.group(1)) - 1;

              // translate for each Ciagr letters a binary code (M=1 and other
              // letters=0)
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

        // get the length of the left outlier
        int leftLengthOutlier = sideWindowsLeft(sequenceCigarBinary);
        informationRead.leftLengthOutlier = leftLengthOutlier;

        // get the length of the right outlier
        int rightLengthOutlier = sideWindowsRight(sequenceCigarBinary);
        informationRead.rightLengthOutlier = rightLengthOutlier;

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

    // test if the sequenceCigarBinary is null
    if (sequenceCigarBinary == null) {
      return 0;
    }

    String window;
    int length = 0;

    // dash of the window in the binary sequence, one by one base
    for (int i = 0; i <= sequenceCigarBinary.length(); i++) {

      // test if the window have finish the sequence
      if (i == (sequenceCigarBinary.length()
          - this.lengthWindowsSideWindow - 2)) {
        break;
      }

      // get the sequence binary on the window
      window =
          sequenceCigarBinary.substring(i, i + this.lengthWindowsSideWindow);

      // test if the sum of this sequence window is superior to the threshold
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

    // test is the binary sequence is null
    if (sequenceCigarBinary == null) {
      return 0;
    }

    String window;
    int length = 0;

    // dash of the window in the binary sequence, one by one base
    for (int i = sequenceCigarBinary.length(); i >= 0; i--) {

      // test if the window have finish the sequence
      if (i == this.lengthWindowsSideWindow) {
        break;
      }

      // get the sequence binary on the window
      window =
          sequenceCigarBinary.substring(i - this.lengthWindowsSideWindow, i);

      // test if the sum of this sequence window is superior to the threshold
      if (sumWindowCIGAR(window) >= this.thresholdSideWindow) {

        // test if the window have finish the sequence
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

    // get each base and add them to make a sum
    for (int i = 0; i <= arrayWindows.length - 1; i++) {

      // test if the base equal to 1
      if (arrayWindows[i] == 1) {
        sum += 1;
      }
    }
    return (double) sum / windows.length();
  }
}
