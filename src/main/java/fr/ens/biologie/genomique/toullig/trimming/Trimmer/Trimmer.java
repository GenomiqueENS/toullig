package fr.ens.biologie.genomique.toullig.trimming.Trimmer;

import java.io.IOException;

/**
 * Interface for Trimmer Created by birer on 03/04/17.
 * @author Aurelien Birer
 */
public interface Trimmer {

  // Method to pre-process trim ONT data
  void preProcessSequence(int leftLengthOutlier, int rightLengthOutlier,
      String sequence, String id, String quality) throws IOException;

  // Method to trim ONT data
  void trim();

}
