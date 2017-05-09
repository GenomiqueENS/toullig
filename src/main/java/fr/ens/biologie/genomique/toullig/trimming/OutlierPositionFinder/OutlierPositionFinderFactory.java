package fr.ens.biologie.genomique.toullig.trimming.OutlierPositionFinder;

import fr.ens.biologie.genomique.toullig.trimming.InformationRead;

import java.io.File;
import java.util.Map;

/**
 * Created by birer on 09/05/17.
 * @author Aurelien Birer
 */
public class OutlierPositionFinderFactory {

  public static OutlierPositionFinder newOutlierPositionFinder(
      Map<String, InformationRead> workTrimmingMap, int addIndexOutlier,
      File fastqFile, int lengthWindowSideWindow, double thresholdSideWindow,
      String mode) {

    OutlierPositionFinder outlierPositionFinder = null;

    // test to process the perfect method
    if (mode.contains("p")) {

      // call PerfectOutlierPositionFinder constructor
      outlierPositionFinder = new PerfectOutlierPositionFinder(workTrimmingMap,
          addIndexOutlier, fastqFile);
    }

    if (mode.contains("sw")) {

      // call SideWindowOutlierPositionFinder constructor
      outlierPositionFinder =
          new SideWindowOutlierPositionFinder(lengthWindowSideWindow,
              thresholdSideWindow, workTrimmingMap, fastqFile);

    }

    if (!mode.contains("sw") && !mode.contains("p")) {

      System.out.println(
          "No Outlier Position Finder method is correctly detected ! Please enter a correct Outlier Position Finder method.");
    }

    return outlierPositionFinder;

  }

}
