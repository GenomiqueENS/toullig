package fr.ens.biologie.genomique.toullig.trimming;

import java.util.Objects;

/**
 * Class of the Working trim Map Created by birer on 03/04/17.
 * @author Aurelien Birer
 */
public class InformationRead {

  public String sequence;
  public String quality;
  public final String cigar;
  public int leftLengthOutlier;
  public int rightLengthOutlier;
  public final int qFlag;
  public final int cigarLength;

  /**
   * Constructor of the InformationRead class.
   * @param sequence, a string sequence of the fastq
   * @param quality, a string quality of the fastq
   * @param cigar, a string cigar of a sam file
   * @param leftLengthOutlier, an int length of the left outlier
   * @param rightLengthOutlier, an int length of the right outlier
   * @param qFlag, an int of the qFlag of a sam file
   * @param cigarLength, an int of the length of the cigar of a sam
   */
  public InformationRead(String sequence, String quality, String cigar,
      int leftLengthOutlier, int rightLengthOutlier, int qFlag,
      int cigarLength) {

    this.sequence = sequence;
    this.quality = quality;
    this.cigar = cigar;
    this.leftLengthOutlier = leftLengthOutlier;
    this.rightLengthOutlier = rightLengthOutlier;
    this.qFlag = qFlag;
    this.cigarLength = cigarLength;

  }

  /**
   * Method of the InformationRead class to get the hashCode.
   * @return , an int of the hashCode
   */
  @Override
  public int hashCode() {
    return Objects.hash(this.sequence, this.quality, this.cigar,
        this.leftLengthOutlier, this.rightLengthOutlier, this.qFlag,
        this.cigarLength);
  }

  /**
   * Method of the InformationRead class to get the fields in string.
   * @return , a string of the InformationRead Object.
   */
  @Override
  public String toString() {
    return this.sequence
        + " " + this.quality + " " + this.cigar + " " + this.leftLengthOutlier
        + " " + rightLengthOutlier + " " + qFlag + " " + cigarLength;
  }

  /**
   * Method of the InformationRead class to compare InformationRead objects.
   * @return , a boolean of the InformationRead Object.
   */
  public boolean equals(InformationRead informationRead) {

    // test if the hashcode is the same between the two InformationRead object
    if (informationRead.hashCode() == Objects.hash(this.sequence, this.quality,
        this.cigar, this.leftLengthOutlier, this.rightLengthOutlier, this.qFlag,
        this.cigarLength)) {

      return true;

    } else {

      return false;

    }

  }

}
