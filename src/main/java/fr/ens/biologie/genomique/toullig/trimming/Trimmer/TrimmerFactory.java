package fr.ens.biologie.genomique.toullig.trimming.Trimmer;

import fr.ens.biologie.genomique.toullig.trimming.InformationRead;

import java.io.File;
import java.io.IOException;
import java.util.Map;

/**
 * Created by birer on 09/05/17.
 * @author Aurelien Birer
 */
public class TrimmerFactory {

  public static Trimmer newTrimmer(Map<String, InformationRead> workTrimmingMap,
      File outputFastqFile, File outputTrimLeftFastaFile,
      File outputTrimRightFastaFile, String adaptorRT,
      String adaptorStrandSwitching, double errorRateCutadapt,
      File fastaLeftOutlierFile, File fastaRightOutlierFile,
      File infoTrimLeftFile, File infoTrimRightFile, File adaptorFile,
      int seedMismatchesTrimmomatic, int palindromeClipThresholdTrimmomatic,
      int simpleClipThreshold, String trimmer, boolean processStatsCutadapt)
      throws IOException {

    Trimmer trimmerObject = null;

    // test to process the cutadapt trimmer
    if (trimmer.contains("cutadapt")) {

      System.out.println("error rate cutadapt: " + errorRateCutadapt);

      // call CutadaptTrimmer constructor
      trimmerObject = new CutadaptTrimmer(workTrimmingMap, outputFastqFile,
          outputTrimLeftFastaFile, outputTrimRightFastaFile, adaptorRT,
          adaptorStrandSwitching, errorRateCutadapt, fastaLeftOutlierFile,
          fastaRightOutlierFile, infoTrimLeftFile, infoTrimRightFile);

      // test to execute stats for cutadapt trim
      if (processStatsCutadapt) {

        // call CutadaptTrimmer constructor
        CutadaptTrimmer trimmingCutadapt =
            new CutadaptTrimmer(workTrimmingMap, outputFastqFile,
                outputTrimLeftFastaFile, outputTrimRightFastaFile, adaptorRT,
                adaptorStrandSwitching, errorRateCutadapt, fastaLeftOutlierFile,
                fastaRightOutlierFile, infoTrimLeftFile, infoTrimRightFile);

        // execute for the left outlier the stats of cutadapt
        trimmingCutadapt.statsLogCutadapt(infoTrimLeftFile,
            "Start stat Left outlier");

        // execute for the right outlier the stats of cutadapt
        trimmingCutadapt.statsLogCutadapt(infoTrimRightFile,
            "Start stat Right outlier");
      }

    }

    // test to process the trimmomatic trimmer
    if (trimmer.contains("trimmomatic")) {

      // call TrimmomaticTrimmer constructor
      trimmerObject =
          new TrimmomaticTrimmer(adaptorFile, seedMismatchesTrimmomatic,
              palindromeClipThresholdTrimmomatic, simpleClipThreshold);
    }

    // test to process no trimmer
    if (trimmer.contains("no")) {

      // call NoTrimmer constructor
      trimmerObject = new NoTrimmer(workTrimmingMap, outputFastqFile);
    }

    return trimmerObject;

  }

}
