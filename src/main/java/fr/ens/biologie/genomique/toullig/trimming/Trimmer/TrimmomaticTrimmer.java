package fr.ens.biologie.genomique.toullig.trimming.Trimmer;

import fr.ens.biologie.genomique.eoulsan.bio.ReadSequence;
import fr.ens.biologie.genomique.eoulsan.bio.io.FastqWriter;
import org.usadellab.trimmomatic.fastq.FastqRecord;
import org.usadellab.trimmomatic.trim.IlluminaClippingTrimmer;
import org.usadellab.trimmomatic.util.Logger;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import static fr.ens.biologie.genomique.toullig.Utils.reverse;
import static fr.ens.biologie.genomique.toullig.trimming.UtilsTrimming.*;

/**
 * Class to execute the Trimmomatic trimmer. Created by birer on 27/03/17.
 * @author Aurelien Birer
 */
public class TrimmomaticTrimmer implements Trimmer {

  private static org.usadellab.trimmomatic.trim.Trimmer trimmer;
  private final File nameOutputFastq;

  public TrimmomaticTrimmer(File adaptorFile, File nameOutputFastq,
      int seedMismatchesTrimmomatic, int palindromeClipThresholdTrimmomatic,
      int simpleClipThreshold) {

    this.nameOutputFastq = nameOutputFastq;

    System.out.println(adaptorFile.getPath()
        + ":" + seedMismatchesTrimmomatic + ":"
        + palindromeClipThresholdTrimmomatic + ":" + simpleClipThreshold);

    try {
      // create the logger for trimmomatic
      Logger logger = new Logger(true, true, true);

      // create the trimmer for trimmomatic
      trimmer = IlluminaClippingTrimmer.makeIlluminaClippingTrimmer(logger,
          adaptorFile.getPath()
              + ":" + seedMismatchesTrimmomatic + ":"
              + palindromeClipThresholdTrimmomatic + ":" + simpleClipThreshold);

    } catch (IOException e) {
      e.printStackTrace();
    }
  }

  //
  // Trimmomatic
  //

  /**
   * Method of the class TrimmomaticTrimmer to trim with the Trimmer interface.
   * @param leftLengthOutlier , the length of the left outlier
   * @param rightLengthOutlier , the length of the left outlier
   * @param sequence, the sequence of the read
   * @param id , the id of the read
   * @param quality , the quality of the read
   */
  public void preProcessSequence(int leftLengthOutlier, int rightLengthOutlier,
      String sequence, String id, String quality) {

    // execute the open of the output fastq file
    try (BufferedWriter fastqOutputFile =
        new BufferedWriter(new FileWriter(this.nameOutputFastq))) {

      //
      // Pre-process before trimmomatic
      //

      // get the left outlier sequence
      String leftOutlierSequence =
          getOutlierLeftSequence(rightLengthOutlier, sequence);

      // get the right outlier sequence
      String rightOutlierSequence =
          getOutlierRightSequence(leftLengthOutlier, sequence);

      // get the left outlier quality
      String leftOutlierQuality =
          getOutlierLeftQuality(rightLengthOutlier, quality);

      // get the right outlier quality
      String rightOutlierQuality =
          getOutlierRightQuality(leftLengthOutlier, quality);

      // get the main sequence (without outlier)
      String mainSequence = sequence.substring(rightLengthOutlier,
          sequence.length() - leftLengthOutlier);

      // get the main quality (without outlier)
      String mainQuality = quality.substring(rightLengthOutlier,
          quality.length() - leftLengthOutlier);

      //
      // Execute trimmomatic
      //

      // Trimming with trimmomatic the right outlier
      String leftTrimSequence = trimmingTrimmomatic(
          reverse(leftOutlierSequence), reverse(leftOutlierQuality));

      // Trimming with trimmomatic the right outlier
      String rigthTrimSequence =
          trimmingTrimmomatic(rightOutlierSequence, rightOutlierQuality);

      //
      // Post-process after trimmomatic
      //

      // get the left trimmed quality
      String leftTrimQuality = quality.substring(
          leftOutlierQuality.length() - leftTrimSequence.length(),
          leftOutlierQuality.length());

      // get the right trimmed quality
      String rigthTrimQuality = quality.substring(mainQuality.length(),
          mainQuality.length() + rightOutlierSequence.length());

      // get the sequence trimmed
      String sequenceTrimmed =
          leftTrimSequence + mainSequence + rigthTrimSequence;

      // get the quality trimmed
      String qualityTrimmed = leftTrimQuality + mainQuality + rigthTrimQuality;

      //
      // Write the sequence and the quality trimmed
      //

      ReadSequence fastq = new ReadSequence();
      fastq.setName(id);
      fastq.setSequence(sequenceTrimmed);
      fastq.setQuality(qualityTrimmed);

      FastqWriter fastqWriter = new FastqWriter(fastqOutputFile);
      fastqWriter.write(fastq);
      fastqWriter.close();

    } catch (IOException e) {
      e.printStackTrace();
    }
  }

  /**
   * Method of the class TrimmomaticTrimmer to trim a sequence with trimmomatic.
   * @return , a sequence trimmed in string format
   */
  String trimmingTrimmomatic(String sequence, String quality) {

    // get basic information of the read
    FastqRecord record = new FastqRecord("name", sequence, "", quality, 33);

    // trim with trimmomatic
    FastqRecord[] result = trimmer.processRecords(new FastqRecord[] {record});

    // return the trimmed sequence
    return result[0].getSequence();
  }

  /**
   * Method of the class TrimmomaticTrimmer to trim with the Trimmer interface
   */
  public void trim() {
  }

}
