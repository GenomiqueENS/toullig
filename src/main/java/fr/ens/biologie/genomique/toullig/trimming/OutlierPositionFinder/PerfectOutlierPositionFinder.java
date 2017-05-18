package fr.ens.biologie.genomique.toullig.trimming.OutlierPositionFinder;

import fr.ens.biologie.genomique.eoulsan.bio.ReadSequence;
import fr.ens.biologie.genomique.eoulsan.bio.io.FastaWriter;
import fr.ens.biologie.genomique.eoulsan.bio.io.FastqReader;
import fr.ens.biologie.genomique.toullig.trimming.InformationRead;
import fr.ens.biologie.genomique.toullig.trimming.Trimmer.Trimmer;

import java.io.File;
import java.io.IOException;
import java.util.Map;

import static fr.ens.biologie.genomique.eoulsan.EoulsanLogger.getLogger;

/**
 * Class to execute the PerfectOutlierPositionFinder class to find outliers
 * Created by birer on 29/03/17.
 * @author Aurelien Birer
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

    // test if the workTrimmingMap Hash is null
    if (workTrimmingMap != null) {

      this.workTrimmingMap = workTrimmingMap;

    } else {

      getLogger().info(
          "Work Trimming Hash is null ! Critical error, contact developpers !");
      this.workTrimmingMap = null;
      System.exit(0);

    }

    this.addIndexOutlier = addIndexOutlier;

    // test if the fastqFile File is null
    if (fastqFile != null) {

      this.fastqFile = fastqFile;

    } else {

      getLogger()
          .info("Fastq File is null ! Critical error, contact developpers !");
      this.fastqFile = null;
      System.exit(0);

    }
  }

  /**
   * Method of the class PerfectOutlierPositionFinder to trim sequence to create
   * sequences files for cutadapt.
   * @param fastaLeftOutlierFile , a fasta File for left outlier
   * @param fastaRightOutlierFile , a fasta File for right outlier
   * @param trimmer, the trimmer interface
   */
  public void findOutliers(File fastaLeftOutlierFile,
      File fastaRightOutlierFile, Trimmer trimmer, File fastqOutputFile)
      throws IOException {

    int countSamReads = 0;
    int countLeftOutlierFind = 0;
    int countRightOutlierFind = 0;
    int countCigarReads = 0;
    int countQFlag16 = 0;
    int countQFlag0 = 0;
    int countQFlag4 = 0;

    // test if the fasta left file is correct
    if (!fastaLeftOutlierFile.isFile() || fastaLeftOutlierFile == null) {

      getLogger().info("Left Fasta File is null or dont exist !");
      System.exit(0);

    }

    // test if the fasta left file is correct
    if (!fastaRightOutlierFile.isFile() || fastaRightOutlierFile == null) {

      getLogger().info("Right Fasta File is null or dont exist !");
      System.exit(0);

    }

    FastaWriter fastqOutputWriter = null;
    FastaWriter fastaLeftOutlierWriter = null;
    FastaWriter fastaRightOutlierWriter = null;

    // test if the trimmer is trimmomatic
    if (trimmer.toString().contains("TrimmomaticTrimmer")) {

      fastqOutputWriter = new FastaWriter(fastqOutputFile);
    }

    // test if the trimmer is cutadapt
    if (trimmer.toString().contains("CutadaptTrimmer")) {

      fastaLeftOutlierWriter = new FastaWriter(fastaLeftOutlierFile);
      fastaRightOutlierWriter = new FastaWriter(fastaRightOutlierFile);
    }

    // open the fastq File
    try (FastqReader reader = new FastqReader(this.fastqFile)) {

      // read the fastq file
      for (ReadSequence read : reader) {

        // get header
        String header = read.getName();

        // get id
        String id = header.substring(0, header.indexOf(" "));

        // get sequence
        String sequence = read.getSequence();

        // get quality
        String quality = read.getQuality();

        // get the information read of the id corresponding in the work trim
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

          // test if the trimmer is trimmomatic
          if (trimmer.toString().contains("TrimmomaticTrimmer")) {

            // pre-process trim
            trimmer.preProcessSequence(leftLengthOutlier, rightLengthOutlier,
                sequence, id, quality, fastqOutputWriter);
          }

          // test if the trimmer is cutadapt
          if (trimmer.toString().contains("CutadaptTrimmer")) {

            // pre-process trim
            trimmer.preProcessSequence(leftLengthOutlier, rightLengthOutlier,
                sequence, id, quality, fastaLeftOutlierWriter,
                fastaRightOutlierWriter);
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

    // test if the trimmer is trimmomatic
    if (trimmer.toString().contains("TrimmomaticTrimmer")) {

      fastqOutputWriter.close();
    }

    // test if the trimmer is cutadapt
    if (trimmer.toString().contains("CutadaptTrimmer")) {

      fastaLeftOutlierWriter.close();
      fastaRightOutlierWriter.close();
    }

    //
    // Display information on the Perfect Outlier Position finder
    //

    getLogger().info("Number of reads in SAM File: " + countSamReads);
    getLogger().info("Number of CIGAR code find (not '*'): " + countCigarReads);
    getLogger().info("Number of QFlag '16': " + countQFlag16);
    getLogger().info("Number of QFlag '0': " + countQFlag0);
    getLogger().info("Number of QFlag '4': " + countQFlag4);
    getLogger().info("Number of left Outlier find: " + countLeftOutlierFind);
    getLogger().info("Number of right Outlier find: " + countRightOutlierFind);

  }
}
