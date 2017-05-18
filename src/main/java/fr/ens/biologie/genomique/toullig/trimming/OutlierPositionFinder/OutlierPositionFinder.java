package fr.ens.biologie.genomique.toullig.trimming.OutlierPositionFinder;

import fr.ens.biologie.genomique.toullig.trimming.Trimmer.Trimmer;

import java.io.File;
import java.io.IOException;

/**
 * Interface for Outlier position finder Created by birer on 31/03/17.
 * @author Aurelien Birer
 */
public interface OutlierPositionFinder {

  /**
   * Method of the interface OutlierPositionFinder to set the Trimmer to use.
   * @param fastaLeftOutlierFile, a FASTA file of the left outlier
   * @param fastaRightOutlierFile, a FASTA file of the right outlier
   * @param trimmer, a trimmer interface to know the trimmer to use
   */
  void findOutliers(File fastaLeftOutlierFile, File fastaRightOutlierFile,
      Trimmer trimmer, File fastqOutputFile) throws IOException;
}
