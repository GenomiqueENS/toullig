package fr.ens.biologie.genomique.toullig.trimming;

import java.io.File;

/**
 * Interface for Outlier position finder Created by birer on 31/03/17.
 */
public interface OutlierPositionFinder {
  // Method to find Outlier in ONT data
  void findOutliers(File fastaLeftOutlierFile, File fastaRightOutlierFile,
      Trimmer trimmer);

}
