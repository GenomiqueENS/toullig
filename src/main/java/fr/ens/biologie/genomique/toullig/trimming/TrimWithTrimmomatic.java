package fr.ens.biologie.genomique.toullig.trimming;

import org.usadellab.trimmomatic.fastq.FastqRecord;
import org.usadellab.trimmomatic.trim.IlluminaClippingTrimmer;
import org.usadellab.trimmomatic.trim.Trimmer;
import org.usadellab.trimmomatic.util.Logger;

import java.io.File;
import java.io.IOException;

/**
 * Created by birer on 27/03/17.
 */
public class TrimWithTrimmomatic {

    private File adaptorFile;
    private int seedMismatchesTrimmomatic;
    private int palindromeClipThresholdTrimmomatic;
    private int simpleClipThreshold;

    public TrimWithTrimmomatic(File adaptorFile, int seedMismatchesTrimmomatic, int palindromeClipThresholdTrimmomatic, int simpleClipThreshold){

        this.adaptorFile=adaptorFile;
        this.seedMismatchesTrimmomatic=seedMismatchesTrimmomatic;
        this.palindromeClipThresholdTrimmomatic=palindromeClipThresholdTrimmomatic;
        this.simpleClipThreshold=simpleClipThreshold;
    }


    //
    // Trimmomatic
    //

    public String TrimmomaticTrim(String sequence, String score) throws IOException {

        System.out.println(this.adaptorFile.getPath()+":"+this.seedMismatchesTrimmomatic+":"+this.palindromeClipThresholdTrimmomatic+":"+this.simpleClipThreshold);
        Logger logger = new Logger(true,true,true);
        Trimmer trimer = IlluminaClippingTrimmer.makeIlluminaClippingTrimmer(logger, this.adaptorFile.getPath()+":"+this.seedMismatchesTrimmomatic+":"+this.palindromeClipThresholdTrimmomatic+":"+this.simpleClipThreshold);
        FastqRecord record = new FastqRecord("name", sequence, "", score, 33);
        FastqRecord [] result = trimer.processRecords(new FastqRecord[] {record});
        return result[0].getSequence();
    }


}
