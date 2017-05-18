package fr.ens.biologie.genomique.toullig.trimming.Trimmer;

import fr.ens.biologie.genomique.eoulsan.bio.ReadSequence;
import fr.ens.biologie.genomique.eoulsan.bio.io.FastaWriter;
import fr.ens.biologie.genomique.eoulsan.bio.io.FastqWriter;
import fr.ens.biologie.genomique.toullig.trimming.InformationRead;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Map;

import static fr.ens.biologie.genomique.eoulsan.EoulsanLogger.getLogger;

/**
 * Class of NoTrimmer method to dont trim sequence
 * @author Aurelien Birer
 */
public class NoTrimmer implements Trimmer {

  private final Map<String, InformationRead> workTrimmingMap;
  private final File nameOutputFastq;

  /**
   * Constructor of the class NoTrimmer
   * @param workTrimmingMap, a Map
   * @param outputFastqFile, a File of output fastq
   */
  public NoTrimmer(Map<String, InformationRead> workTrimmingMap,
      File outputFastqFile) {

    this.workTrimmingMap = workTrimmingMap;
    this.nameOutputFastq = outputFastqFile;
  }

  /**
   * Method of the class NoTrimmer to write the sequence
   */
  private void writeSequenceWithNoOutlier() {

    int i = 0;
    int countWritten = 0;
    int countNull = 0;
    String shortestFastqSequence = "";

    // test if the output left fasta file is correctly open
    try (BufferedWriter fastqTrimFile =
        new BufferedWriter(new FileWriter(this.nameOutputFastq))) {

      // write the trimmed read
      FastqWriter fastqWriter = new FastqWriter(fastqTrimFile);

      // get id for each red on the work map
      for (String id : this.workTrimmingMap.keySet()) {

        i++;

        // get information for the read
        InformationRead informationRead = this.workTrimmingMap.get(id);

        String cigar = informationRead.cigar;
        String quality = informationRead.quality;
        String sequence = informationRead.sequence;
        int leftLengthOutlier = informationRead.leftLengthOutlier;
        int rightLengthOutlier = informationRead.rightLengthOutlier;

        // test if the read is unmapped
        if (!cigar.equals("*")) {

          // test if the length of the left outlier is negatif
          if (leftLengthOutlier < 0) {
            leftLengthOutlier = 0;
          }

          // test if the length of the left outlier is negatif
          if (rightLengthOutlier < 0) {
            rightLengthOutlier = 0;
          }

          String mainSequenceWithoutOutlier = "";

          String trimmedSequence = "";
          String trimmedQuality = "";

          // test if the rightlengthsequence and the leftlengthsequence are
          // overlap
          if ((sequence.length() - rightLengthOutlier) > leftLengthOutlier) {
            trimmedSequence = sequence.substring(leftLengthOutlier,
                sequence.length() - rightLengthOutlier);
            trimmedQuality = quality.substring(leftLengthOutlier,
                quality.length() - rightLengthOutlier);
          }

          // test if the quality length is differerent to the sequence length
          if (trimmedQuality.length() != trimmedSequence.length()) {
            getLogger().info("problem :  "
                + trimmedQuality.length() + "     " + trimmedSequence.length()
                + "   " + mainSequenceWithoutOutlier.length());
            getLogger().info(leftLengthOutlier
                + "     " + rightLengthOutlier + "     " + id);
          }

          // test if the sequence trimmed is empty
          if (!trimmedSequence.isEmpty()) {

            ReadSequence fastq = new ReadSequence();
            fastq.setName(id);
            fastq.setSequence(trimmedSequence);
            fastq.setQuality(trimmedQuality);

            // writer the fastq read
            fastqWriter.write(fastq);

            countWritten++;

            // test for the first loop
            if (i == 1) {
              shortestFastqSequence = trimmedSequence;
            }

            // test if a new shortest read is comput in this loop
            if (shortestFastqSequence.length() >= trimmedSequence.length()) {
              shortestFastqSequence = trimmedSequence;
            }

          } else {
            countNull++;
          }
        }
      }

      // close the writer
      fastqTrimFile.close();

    } catch (IOException e) {
      e.printStackTrace();
    }

    getLogger()
        .info("The shortest read size is: " + shortestFastqSequence.length());
    getLogger().info("Number of trim read write: " + countWritten);
    getLogger().info("Number of trim read null: " + countNull);

  }

  // No preprocessTrimming are require
  @Override
  public void preProcessSequence(int leftLengthOutlier, int rightLengthOutlier,
      String sequence, String id, String quality,
      FastaWriter fastaLeftOutlierWriter, FastaWriter fastaRightOutlierWriter) {

  }

  @Override
  public void preProcessSequence(int leftLengthOutlier, int rightLengthOutlier,
      String sequence, String id, String quality, FastaWriter fastqOutputFile)
      throws IOException {

  }

  /**
   * Method of the class NoTrimmer to trim with the Trimmer interface
   */
  public void trim() {

    // execute the writting of the sequence.
    writeSequenceWithNoOutlier();

  }

}