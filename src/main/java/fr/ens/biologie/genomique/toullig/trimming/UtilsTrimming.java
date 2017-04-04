package fr.ens.biologie.genomique.toullig.trimming;

import java.io.IOException;

import fr.ens.biologie.genomique.eoulsan.bio.Sequence;
import fr.ens.biologie.genomique.eoulsan.bio.io.FastaWriter;

/**
 * Class of somes utils methods for cutadaptTrimming. Created by birer on
 * 29/03/17.
 */
class UtilsTrimming {

  /**
   * Method of the class TrimModes to obtain the sequence of the left outlier.
   * @param lengthOutlierBegin, the length of the outlier
   * @param sequence, the sequence of the read
   * @return the sequence of the left outlier
   */
  static String getOutlierLeftSequence(int lengthOutlierBegin,
      String sequence) {
    return sequence.substring(0, lengthOutlierBegin);
  }

  /**
   * Method of the class TrimModes to obtain the sequence of the right outlier.
   * @param lengthOutlierEnd, the length of the outlier
   * @param sequence, the sequence of the read
   * @return the sequence of the right outlier
   */
  static String getOutlierRightSequence(int lengthOutlierEnd, String sequence) {

    return sequence.substring(sequence.length() - lengthOutlierEnd,
        sequence.length());
  }

  /**
   * Method of the class TrimModes to obtain the score of the left outlier.
   * @param lengthOutlierBegin, the length of the outlier
   * @param score, the score of the read
   * @return , a string of quality
   */
  static String getOutlierLeftQuality(int lengthOutlierBegin, String score) {
    return score.substring(0, lengthOutlierBegin);
  }

  /**
   * Method of the class TrimModes to obtain the score of the right outlier.
   * @param lengthOutlierEnd, the length of the outlier
   * @param score, the score of the read
   * @return , a string of quality
   */
  static String getOutlierRightQuality(int lengthOutlierEnd, String score) {
    return score.substring(score.length() - lengthOutlierEnd, score.length());
  }

  //
  // get Outlier
  //

  /**
   * Method of the class TrimModes to write the outliers (3' and 5') in fasta
   * files.
   * @param leftLengthOutlier, the length of the outlier right
   * @param rightLengthOutlier, the length of the outlier left
   * @param sequence, the sequence of the read
   * @param id, the id of the read
   */
  static void writeOutliers(int leftLengthOutlier, int rightLengthOutlier,
      String sequence, String id, FastaWriter leftFastaWriter,
      FastaWriter rightFastaWriter) {

    // get the left outlier
    String leftOutlierSequence =
        getOutlierLeftSequence(leftLengthOutlier, sequence);

    // get the right outlier
    String rightOutlierSequence =
        getOutlierRightSequence(rightLengthOutlier, sequence);

    try {

      Sequence leftFasta = new Sequence();

      // get id and sequence of the left outlier
      leftFasta.setName(id);
      leftFasta.setSequence(leftOutlierSequence);

      // write the left outlier
      leftFastaWriter.write(leftFasta);

      Sequence rightFasta = new Sequence();

      // get id and sequence of the right outlier
      rightFasta.setName(id);
      rightFasta.setSequence(rightOutlierSequence);

      // write the right outlier
      rightFastaWriter.write(rightFasta);

    } catch (IOException e) {
      e.printStackTrace();
    }
  }
}
