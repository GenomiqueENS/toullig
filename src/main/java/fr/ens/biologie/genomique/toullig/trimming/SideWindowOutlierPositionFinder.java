package fr.ens.biologie.genomique.toullig.trimming;

import fr.ens.biologie.genomique.eoulsan.bio.ReadSequence;
import fr.ens.biologie.genomique.eoulsan.bio.io.FastaWriter;
import fr.ens.biologie.genomique.eoulsan.bio.io.FastqWriter;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static fr.ens.biologie.genomique.toullig.Utils.*;
import static fr.ens.biologie.genomique.toullig.trimming.UtilsTrimming.*;

/**
 * Class to execute the SideWindow method to find outliers. Created by birer on
 * 29/03/17.
 */
public class SideWindowOutlierPositionFinder {

  private int lengthWindowsSideWindow;
  private double thresholdSideWindow;
  private Map<String, String[]> workTrimmingMap;
  private File nameOutputFastq;
  private File adaptorFile;
  private int seedMismatchesTrimmomatic;
  private int palindromeClipThresholdTrimmomatic;
  private int simpleClipThreshold;
  private boolean processCutadapt;
  private boolean processTrimmomatic;

  public SideWindowOutlierPositionFinder(int lengthWindowsSideWindow,
      double thresholdSideWindow, Map<String, String[]> workTrimmingMap,
      File nameOutputFastq, File adaptorFile, int seedMismatchesTrimmomatic,
      int palindromeClipThresholdTrimmomatic, int simpleClipThreshold,
      boolean processCutadapt, boolean processTrimmomatic) throws IOException {

    this.lengthWindowsSideWindow = lengthWindowsSideWindow;
    this.thresholdSideWindow = thresholdSideWindow;
    this.workTrimmingMap = workTrimmingMap;
    this.nameOutputFastq = nameOutputFastq;
    this.adaptorFile = adaptorFile;
    this.seedMismatchesTrimmomatic = seedMismatchesTrimmomatic;
    this.palindromeClipThresholdTrimmomatic =
        palindromeClipThresholdTrimmomatic;
    this.simpleClipThreshold = simpleClipThreshold;
    this.processCutadapt = processCutadapt;
    this.processTrimmomatic = processTrimmomatic;
  }

  /**
   * Method of the class SideWindowOutlierPositionFinder to trim sequence to
   * create sequences files for cutadapt with the side-window method.
   * @throws IOException
   * @throws InterruptedException
   */
  public void getOutliers(File fastaLeftOutlierFile, File fastaRightOutlierFile)
      throws IOException, InterruptedException {

    Pattern oneCodeCigarPattern = Pattern.compile("(([0-9]*[A-Z]).*)");
    Pattern lengthOneCodeCigarPattern = Pattern.compile("([0-9]*)(.)");
    int lengthSequenceCigar;

    try (
        BufferedWriter fastqOutputFile =
            new BufferedWriter(new FileWriter(this.nameOutputFastq));
        FastaWriter leftFastaWriter = new FastaWriter(fastaLeftOutlierFile);
        FastaWriter rightFastaWriter = new FastaWriter(fastaRightOutlierFile)) {

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

          if (this.processCutadapt) {
            writeOutliers(leftLengthOutlier, rightLengthOutlier, sequence, id,
                leftFastaWriter, rightFastaWriter);
          }

          if (this.processTrimmomatic) {

            String leftOutlierSequence =
                getOutlierLeftSequence(leftLengthOutlier, sequence);
            String rightOutlierSequence =
                getOutlierRightSequence(rightLengthOutlier, sequence);
            String leftOutlierQuality =
                getOutlierLeftQuality(leftLengthOutlier, quality);
            String rightOutlierQuality =
                getOutlierRightQuality(rightLengthOutlier, quality);

            String mainSequence = sequence.substring(leftLengthOutlier,
                sequence.length() - rightLengthOutlier);
            String mainQuality = quality.substring(leftLengthOutlier,
                quality.length() - rightLengthOutlier);

            TrimmomaticTrimmer trimmingTrimmomatic = new TrimmomaticTrimmer(
                this.adaptorFile, this.seedMismatchesTrimmomatic,
                this.palindromeClipThresholdTrimmomatic,
                this.simpleClipThreshold);

            String rigthTrimSequence = trimmingTrimmomatic
                .TrimmomaticTrim(rightOutlierSequence, rightOutlierQuality);
            String leftTrimSequence = trimmingTrimmomatic.TrimmomaticTrim(
                reverse(leftOutlierSequence), reverse(leftOutlierQuality));

            String rigthTrimQuality = quality.substring(mainQuality.length(),
                mainQuality.length() + rightOutlierSequence.length());
            String leftTrimQuality = quality.substring(
                leftOutlierQuality.length() - leftTrimSequence.length(),
                leftOutlierQuality.length());

            String sequenceTranscript =
                leftTrimSequence + mainSequence + rigthTrimSequence;
            String qualityTranscript =
                leftTrimQuality + mainQuality + rigthTrimQuality;

            ReadSequence fastq = new ReadSequence();
            fastq.setName(id);
            fastq.setSequence(sequenceTranscript);
            fastq.setQuality(qualityTranscript);
            FastqWriter fastqWriter = new FastqWriter(fastqOutputFile);
            fastqWriter.write(fastq);
          }
        }
      }
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
