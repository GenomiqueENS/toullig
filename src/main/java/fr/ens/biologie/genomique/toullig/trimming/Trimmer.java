package fr.ens.biologie.genomique.toullig.trimming;

/**
 * Interface for Trimmer Created by birer on 03/04/17.
 */
public interface Trimmer {

  // Method to pre-process trimming ONT data
  void preProcessTrimming(int leftLengthOutlier, int rightLengthOutlier,
      String sequence, String id, String quality);

  // Method to trimming ONT data
  void trimming();

}
