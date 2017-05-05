package fr.ens.biologie.genomique.toullig.trimming.OutlierPositionFinder;

import fr.ens.biologie.genomique.eoulsan.bio.ReadSequence;
import fr.ens.biologie.genomique.eoulsan.bio.io.FastqReader;
import fr.ens.biologie.genomique.toullig.trimming.InformationRead;
import fr.ens.biologie.genomique.toullig.trimming.Trimmer.Trimmer;

import java.io.File;
import java.io.IOException;
import java.util.Map;

/**
 * Class to execute the PerfectOutlierPositionFinder class to find outliers
 * Created by birer on 29/03/17.
 */
public class PerfectOutlierPositionFinder implements OutlierPositionFinder {

  private final Map<String, InformationRead> workTrimmingMap;
  private final int addIndexOutlier;
  private final File fastqFile;

  /**
   * Constructor of the PerfectOutlierPositionFinder class.
   * @param workTrimmingMap, a Map of work information
   * @param addIndexOutlier, a int of outlier index add
   */
  public PerfectOutlierPositionFinder(
      Map<String, InformationRead> workTrimmingMap, int addIndexOutlier,
      File fastqFile) {

    this.workTrimmingMap = workTrimmingMap;
    this.addIndexOutlier = addIndexOutlier;
    this.fastqFile = fastqFile;
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

        int leftLengthOutlier = 0;
        int rightLengthOutlier = 0;
        countSamReads++;

        // get information of each read
        String cigar = informationRead.cigar;
        int qFlag = informationRead.qFlag;

        // trim by CIGAR
        if (!"*".equals(cigar)) {

          // Locate index between the adaptor and the mRNA with the Soft and
          // Hard clipping for Left side
          int leftIndexOutlier = cigar.indexOf("S");

          // test if the last index is a Hard clipping
          if (leftIndexOutlier > cigar.indexOf("H") && cigar.contains("H")) {
            leftIndexOutlier = cigar.indexOf("H");
          }

          // test if the left index outlier is correct
          if (leftIndexOutlier != -1
              && leftIndexOutlier + 1 != cigar.length()
              && leftIndexOutlier <= 6) {

            // get the left length of the outlier
            leftLengthOutlier =
                Integer.parseInt(cigar.substring(0, leftIndexOutlier))
                    + this.addIndexOutlier;

            // strand case
            if (qFlag == 0) {

              informationRead.leftLengthOutlier = leftLengthOutlier;

            }
            // reverse complement case
            else {

              informationRead.rightLengthOutlier = leftLengthOutlier;
            }

            countLeftOutlierFind++;
          } else {
            rightLengthOutlier = 0;
          }

          // Locate index between the adaptor and the mRNA with the Soft and
          // Hard clipping for Right side
          int rightIndexOutlier = cigar.lastIndexOf("S");

          // test if the last index is a Hard clipping
          if (rightIndexOutlier < cigar.lastIndexOf("H")) {
            rightIndexOutlier = cigar.lastIndexOf("H");
          }

          // test if the right index is negatif
          if (rightIndexOutlier == -1) {
            rightIndexOutlier = 0;
          }

          // test if the left index outlier is correct
          if (rightIndexOutlier + 1 == cigar.length()
              && cigar.length() - rightIndexOutlier < 5) {

            // get the right length of the outlier
            rightLengthOutlier = Integer.parseInt(
                cigar.substring(cigar.lastIndexOf("M") + 1, rightIndexOutlier))
                + this.addIndexOutlier;

            // strand case
            if (qFlag == 0) {

              informationRead.rightLengthOutlier = rightLengthOutlier;

            }
            // reverse complement case
            else {

              informationRead.leftLengthOutlier = rightLengthOutlier;
            }

            countRightOutlierFind++;

          } else {
            leftLengthOutlier = 0;
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

        } else {
          countCigarReads++;
        }

        // test if the qFlag equal to 16
        if (qFlag == 16) {
          countQFlag16++;
        }

        // test if the qFlag equal to 0
        if (qFlag == 0) {
          countQFlag0++;
        }

        // test if the qFlag equal to 4
        if (qFlag == 4) {
          countQFlag4++;
        }
      }

      reader.close();
    } catch (IOException e) {
      e.printStackTrace();
    }

    //
    // Display information on the Perfect Outlier Position finder
    //

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
