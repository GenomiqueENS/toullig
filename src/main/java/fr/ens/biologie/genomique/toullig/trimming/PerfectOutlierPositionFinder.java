package fr.ens.biologie.genomique.toullig.trimming;

import java.io.File;
import java.util.Map;

/**
 * Class to execute the PerfectOutlierPositionFinder class to find outliers
 * Created by birer on 29/03/17.
 */
public class PerfectOutlierPositionFinder implements OutlierPositionFinder {

  private final Map<String, String[]> workTrimmingMap;
  private final int addIndexOutlier;

  public PerfectOutlierPositionFinder(Map<String, String[]> workTrimmingMap,
      int addIndexOutlier) {

    this.workTrimmingMap = workTrimmingMap;
    this.addIndexOutlier = addIndexOutlier;
  }

  /**
   * Method of the class PerfectOutlierPositionFinder to trim sequence to create
   * sequences files for cutadapt.
   * @param fastaLeftOutlierFile , a fasta File for left outlier
   * @param fastaRightOutlierFile , a fasta File for right outlier
   * @param trimmer, the trimmer interface
   */
  public void findOutliers(File fastaLeftOutlierFile,
      File fastaRightOutlierFile, Trimmer trimmer) {

    int countSamReads = 0;
    int countLeftOutlierFind = 0;
    int countRightOutlierFind = 0;
    int countCigarReads = 0;
    int countQFlag16 = 0;
    int countQFlag0 = 0;
    int countQFlag4 = 0;

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

        // pre-process trimming
        trimmer.preProcessTrimming(leftLengthOutlier, rightLengthOutlier,
            sequence, id, quality);

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
}
