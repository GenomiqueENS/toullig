package fr.ens.biologie.genomique.toullig.trimming;

import java.io.BufferedWriter;
import java.io.IOException;

import fr.ens.biologie.genomique.toullig.Utils;

/**
 * Created by birer on 29/03/17.
 */
public class UtilsTrimming {

    public UtilsTrimming(){

    }

    /**
     * Method of the class TrimModes to obtain the sequence of the left outlier.
     * @param lengthOutlierBegin, the length of the outlier
     * @param sequence, the sequence of the read
     * @return the sequence of the left outlier
     */
    public String getOutlierLeftSequence(int lengthOutlierBegin, String sequence){
        return sequence.substring(0,lengthOutlierBegin);
    }

    /**
     * Method of the class TrimModes to obtain the sequence of the right outlier.
     * @param lengthOutlierEnd, the length of the outlier
     * @param sequence, the sequence of the read
     * @return the sequence of the right outlier
     */
    public String getOutlierRightSequence(int lengthOutlierEnd, String sequence){

        return sequence.substring(sequence.length()-lengthOutlierEnd,sequence.length());
    }

    /**
     * Method of the class TrimModes to obtain the score of the left outlier.
     * @param lengthOutlierBegin, the length of the outlier
     * @param score, the score of the read
     * @return
     */
    public String getOutlierLeftScore(int lengthOutlierBegin, String score){
        return score.substring(0,lengthOutlierBegin);
    }

    /**
     * Method of the class TrimModes to obtain the score of the right outlier.
     * @param lengthOutlierEnd, the length of the outlier
     * @param score, the score of the read
     * @return
     */
    public String getOutlierRightScore(int lengthOutlierEnd, String score){
        return score.substring(score.length()-lengthOutlierEnd,score.length());
    }

    //
    // get Outlier
    //

    /**
     * Method of the class TrimModes to write the outliers (3' and 5') in fasta files.
     * @param lengthOutlierBegin, the length of the outlier right
     * @param lengthOutlierEnd, the length of the outlier left
     * @param sequence, the sequence of the read
     * @param ID, the ID of the read
     * @param utils, a utils object
     * @throws IOException
     */
    public void writeOutliers(int lengthOutlierBegin, int lengthOutlierEnd, String sequence, String ID, Utils utils, BufferedWriter fastaFileLeftOutlier, BufferedWriter fastaFileRightOutlier) throws IOException {

        String leftOutlierSequence= getOutlierLeftSequence(lengthOutlierBegin,sequence);
        String rightOutlierSequence= getOutlierRightSequence(lengthOutlierEnd,sequence);
        utils.writeFasta(leftOutlierSequence,ID,fastaFileLeftOutlier);
        utils.writeFasta(rightOutlierSequence,ID,fastaFileRightOutlier);
    }
}
