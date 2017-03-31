package fr.ens.biologie.genomique.toullig.trimming;

import fr.ens.biologie.genomique.eoulsan.bio.ReadSequence;
import fr.ens.biologie.genomique.eoulsan.bio.io.FastaWriter;
import fr.ens.biologie.genomique.eoulsan.bio.io.FastqWriter;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Map;

import static fr.ens.biologie.genomique.toullig.Utils.*;
import static fr.ens.biologie.genomique.toullig.trimming.UtilsTrimming.*;

/**
 * Class to execute the PerfectOutlierPositionFinder class to find outliers
 * Created by birer on 29/03/17.
 */
public class PerfectOutlierPositionFinder implements AutoCloseable {

  private Map<String, String[]> workTrimmingMap;
  private File nameOutputFastq;
  private int addIndexOutlier;
  private File adaptorFile;
  private int seedMismatchesTrimmomatic;
  private int palindromeClipThresholdTrimmomatic;
  private int simpleClipThreshold;
  private boolean processCutadapt;
  private boolean processTrimmomatic;

  public PerfectOutlierPositionFinder(Map<String, String[]> workTrimmingMap,
      File nameOutputFastq, int addIndexOutlier, File adaptorFile,
      int seedMismatchesTrimmomatic, int palindromeClipThresholdTrimmomatic,
      int simpleClipThreshold, boolean processCutadapt,
      boolean processTrimmomatic) throws IOException {

    this.workTrimmingMap = workTrimmingMap;
    this.nameOutputFastq = nameOutputFastq;
    this.addIndexOutlier = addIndexOutlier;
    this.adaptorFile = adaptorFile;
    this.seedMismatchesTrimmomatic = seedMismatchesTrimmomatic;
    this.palindromeClipThresholdTrimmomatic =
        palindromeClipThresholdTrimmomatic;
    this.simpleClipThreshold = simpleClipThreshold;
    this.processCutadapt = processCutadapt;
    this.processTrimmomatic = processTrimmomatic;
  }

  /**
   * Method of the class PerfectOutlierPositionFinder to trim sequence to create
   * sequences files for cutadapt.
   */
  public void findOutliers(File fastaLeftOutlierFile,
      File fastaRightOutlierFile) {

    int countSamReads = 0;
    int countLeftOutlierFind = 0;
    int countRightOutlierFind = 0;
    int countCigarReads = 0;
    int countQFlag16 = 0;
    int countQFlag0 = 0;
    int countQFlag4 = 0;

    try (
        BufferedWriter fastqOutputFile =
            new BufferedWriter(new FileWriter(this.nameOutputFastq));
        FastaWriter leftFastaWriter = new FastaWriter(fastaLeftOutlierFile);
        FastaWriter rightFastaWriter = new FastaWriter(fastaRightOutlierFile)) {

      for (String id : this.workTrimmingMap.keySet()) {

        int leftLengthOutlier = 0;
        int rightLengthOutlier = 0;
        countSamReads++;
        String[] tabValue = this.workTrimmingMap.get(id);
        String sequence = tabValue[0];
        String quality = tabValue[1];
        String cigar = tabValue[2];
        int qFlag = Integer.parseInt(tabValue[5]);

        // trim by CIGAR
        if (!"*".equals(cigar)) {

          // Locate index between the adaptor and the mRNA with the Soft and
          // Hard clipping for Left side
          int leftIndexOutlier = cigar.indexOf("S");
          if (leftIndexOutlier > cigar.indexOf("H") && cigar.contains("H")) {
            leftIndexOutlier = cigar.indexOf("H");
          }

          if (leftIndexOutlier != -1
              && leftIndexOutlier + 1 != cigar.length()
              && leftIndexOutlier <= 6) {

            leftLengthOutlier =
                Integer.parseInt(cigar.substring(0, leftIndexOutlier))
                    + this.addIndexOutlier;
            tabValue[3] = "" + leftLengthOutlier;

            countLeftOutlierFind++;
          } else {
            rightLengthOutlier = 0;
          }

          // Locate index between the adaptor and the mRNA with the Soft and
          // Hard clipping for Right side
          int rightIndexOutlier = cigar.lastIndexOf("S");
          if (rightIndexOutlier < cigar.lastIndexOf("H")) {
            rightIndexOutlier = cigar.lastIndexOf("H");
          }
          if (rightIndexOutlier == -1) {
            rightIndexOutlier = 0;
          }

          if (rightIndexOutlier + 1 == cigar.length()
              && cigar.length() - rightIndexOutlier < 5) {

            rightLengthOutlier = Integer.parseInt(
                cigar.substring(cigar.lastIndexOf("M") + 1, rightIndexOutlier))
                + this.addIndexOutlier;
            tabValue[4] = "" + rightLengthOutlier;
            countRightOutlierFind++;
          } else {
            leftLengthOutlier = 0;
          }

          if (this.processCutadapt) {
            writeOutliers(leftLengthOutlier, rightLengthOutlier, sequence, id,
                leftFastaWriter, rightFastaWriter);
          }

          if (this.processTrimmomatic) {

            String leftOutlierSequence =
                getOutlierLeftSequence(rightLengthOutlier, sequence);
            String rightOutlierSequence =
                getOutlierRightSequence(leftLengthOutlier, sequence);
            String leftOutlierScore =
                getOutlierLeftQuality(rightLengthOutlier, quality);
            String rightOutlierScore =
                getOutlierRightQuality(leftLengthOutlier, quality);

            String mainSequence = sequence.substring(rightLengthOutlier,
                sequence.length() - leftLengthOutlier);
            String mainScore = quality.substring(rightLengthOutlier,
                quality.length() - leftLengthOutlier);

            TrimmomaticTrimmer trimmingTrimmomatic = new TrimmomaticTrimmer(
                this.adaptorFile, this.seedMismatchesTrimmomatic,
                this.palindromeClipThresholdTrimmomatic,
                this.simpleClipThreshold);

            String rigthTrimSequence = trimmingTrimmomatic
                .TrimmomaticTrim(rightOutlierSequence, rightOutlierScore);
            String leftTrimSequence = trimmingTrimmomatic.TrimmomaticTrim(
                reverse(leftOutlierSequence), reverse(leftOutlierScore));

            String rigthTrimScore = quality.substring(mainScore.length(),
                mainScore.length() + rightOutlierSequence.length());
            String leftTrimScore = quality.substring(
                leftOutlierScore.length() - leftTrimSequence.length(),
                leftOutlierScore.length());

            String sequenceTranscript =
                leftTrimSequence + mainSequence + rigthTrimSequence;
            String scoreTranscript = leftTrimScore + mainScore + rigthTrimScore;

            ReadSequence fastq = new ReadSequence();
            fastq.setName(id);
            fastq.setSequence(sequenceTranscript);
            fastq.setQuality(scoreTranscript);
            FastqWriter fastqWriter = new FastqWriter(fastqOutputFile);
            fastqWriter.write(fastq);

          }
        } else {
          countCigarReads++;
        }
        if (qFlag == 16) {
          countQFlag16++;
        }
        if (qFlag == 0) {
          countQFlag0++;
        }
        if (qFlag == 4) {
          countQFlag4++;
        }
      }
    } catch (IOException e) {
      e.printStackTrace();
    }

    System.out.println("Number of reads in SAM File: " + countSamReads);
    System.out
        .println("Number of CIGAR code find (not '*'): " + countCigarReads);
    System.out.println("Number of QFlag '16': " + countQFlag16);
    System.out.println("Number of QFlag '0': " + countQFlag0);
    System.out.println("Number of QFlag '4': " + countQFlag4);
    System.out.println("Number of left Outlier find: " + countLeftOutlierFind);
    System.out
        .println("Number of right Outlier find: " + countRightOutlierFind);

  }

  @Override
  public void close() throws Exception {

  }
}
