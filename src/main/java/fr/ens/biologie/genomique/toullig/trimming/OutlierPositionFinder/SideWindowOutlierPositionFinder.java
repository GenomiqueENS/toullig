package fr.ens.biologie.genomique.toullig.trimming.OutlierPositionFinder;

import fr.ens.biologie.genomique.eoulsan.bio.ReadSequence;
import fr.ens.biologie.genomique.eoulsan.bio.io.FastqReader;
import fr.ens.biologie.genomique.toullig.trimming.InformationRead;
import fr.ens.biologie.genomique.toullig.trimming.Trimmer.Trimmer;

import java.io.File;
import java.io.IOException;
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
  private final File fastqFile;

  /**
   * Constructor of the SideWindow class.
   * @param lengthWindowsSideWindow, a int of the length of the window for the
   *          side-window outlier postion finder
   * @param thresholdSideWindow, a double of the threshold for the side-window
   *          outlier postion finder
   * @param workTrimmingMap, a Map of working information
   */
  public SideWindowOutlierPositionFinder(int lengthWindowsSideWindow,
      double thresholdSideWindow, Map<String, InformationRead> workTrimmingMap,
      File fastqFile) {

    this.lengthWindowsSideWindow = lengthWindowsSideWindow;
    this.thresholdSideWindow = thresholdSideWindow;
    this.workTrimmingMap = workTrimmingMap;
    this.fastqFile = fastqFile;
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

    // open the fastq File
    try (FastqReader reader = new FastqReader(this.fastqFile)) {

      // read the fastq file
      for (ReadSequence read : reader) {

        // get header
        String header = read.getName();

        // split header for get id
        String[] part = header.split(" ");

        // get id
        String id = part[0];

        // get sequence
        String sequence = read.getSequence();

        // get quality
        String quality = read.getQuality();

        // get the information read of the id corresponding in the work trimming
        // map
        InformationRead informationRead = this.workTrimmingMap.get(id);

        // case that the SAM file has not succesful read this id (present in
        // fastq file) read
        if (informationRead == null) {

          break;
        }

        // set the sequence of the read
        informationRead.sequence = sequence;

        // set the quality of the read
        informationRead.quality = quality;

        StringBuilder sequenceCigarBinary = new StringBuilder();

        // get information of the read
        String cigar = informationRead.cigar;

        // get the qflag
        int qFlag = informationRead.qFlag;

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
                    sequenceCigarBinary.append(1);
                    continue;
                  }
                  if (lengthOneCodeCigarMatcher.group(2).equals("N")
                      || lengthOneCodeCigarMatcher.group(2).equals("D")) {
                    continue;
                  }
                  if (!lengthOneCodeCigarMatcher.group(2).equals("N")
                      || !lengthOneCodeCigarMatcher.group(2).equals("M")) {
                    sequenceCigarBinary.append(0);
                  }
                }
              }
            }
          }

          int leftLengthOutlier = 0;
          int rightLengthOutlier = 0;

          // strand case
          if (qFlag == 0) {

            // get the length of the left outlier
            leftLengthOutlier = sideWindowsLeft(sequenceCigarBinary.toString());
            informationRead.leftLengthOutlier = leftLengthOutlier;

            // get the length of the right outlier
            rightLengthOutlier =
                sideWindowsRight(sequenceCigarBinary.toString());
            informationRead.rightLengthOutlier = rightLengthOutlier;

            // System.out.println(leftLengthOutlier
            // + " " + rightLengthOutlier + " " + sequence.length());

          }
          // reverse complement case
          else {

            // get the length of the left outlier
            leftLengthOutlier = sideWindowsLeft(sequenceCigarBinary.toString());
            informationRead.rightLengthOutlier = leftLengthOutlier;

            // get the length of the right outlier
            rightLengthOutlier =
                sideWindowsRight(sequenceCigarBinary.toString());
            informationRead.leftLengthOutlier = rightLengthOutlier;

            // System.out.println(leftLengthOutlier
            // + " " + rightLengthOutlier + " " + sequence.length());
          }

          // strand case
          if (qFlag == 0) {

            // pre-process trimming
            trimmer.preProcessTrimming(leftLengthOutlier, rightLengthOutlier,
                sequence, id, quality);

          }
          // reverse complement case
          else {

            // pre-process trimming
            trimmer.preProcessTrimming(rightLengthOutlier, leftLengthOutlier,
                sequence, id, quality);
          }
        }
      }
      reader.close();
    } catch (IOException e) {
      e.printStackTrace();
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
