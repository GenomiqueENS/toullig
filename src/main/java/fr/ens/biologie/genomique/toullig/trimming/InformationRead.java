package fr.ens.biologie.genomique.toullig.trimming;

import java.util.Objects;

/**
 * Class of the Working trimming Map Created by birer on 03/04/17.
 */
public class InformationRead {

  public String sequence;
  public String quality;
  public final String cigar;
  public int leftLengthOutlier;
  public int rightLengthOutlier;
  public final int qFlag;
  public final int cigarLength;

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

  @Override
  public int hashCode() {
    return Objects.hash(this.sequence, this.quality, this.cigar,
        this.leftLengthOutlier, this.rightLengthOutlier, this.qFlag,
        this.cigarLength);
  }

  @Override
  public String toString() {
    return this.sequence
        + " " + this.quality + " " + this.cigar + " " + this.leftLengthOutlier
        + " " + rightLengthOutlier + " " + qFlag + " " + cigarLength;
  }

}
