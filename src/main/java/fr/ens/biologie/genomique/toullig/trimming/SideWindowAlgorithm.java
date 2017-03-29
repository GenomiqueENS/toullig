package fr.ens.biologie.genomique.toullig.trimming;

import fr.ens.biologie.genomique.toullig.Utils;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Created by birer on 29/03/17.
 */
public class SideWindowAlgorithm {

    private int lengthWindowsSW;
    private double thresholdSW;
    private HashMap<String, String[]> fastqHash;
    private File nameOutputFastq;
    private File adaptorFile;
    private int seedMismatchesTrimmomatic;
    private int palindromeClipThresholdTrimmomatic;
    private int simpleClipThreshold;
    private boolean processCutadapt;
    private boolean processTrimmomatic;
    private BufferedWriter fastaFileLeftOutlier;
    private BufferedWriter fastaFileRightOutlier;

    public SideWindowAlgorithm(int lengthWindowsSW, double thresholdSW, HashMap fastqHash, File nameOutputFastq, File adaptorFile, int seedMismatchesTrimmomatic, int palindromeClipThresholdTrimmomatic, int simpleClipThreshold, boolean processCutadapt, boolean processTrimmomatic, String pathFastaFileLeftOutlier, String pathFastaFileRightOutlier) throws IOException {

        this.lengthWindowsSW=lengthWindowsSW;
        this.thresholdSW=thresholdSW;
        this.fastqHash=fastqHash;
        this.nameOutputFastq=nameOutputFastq;
        this.adaptorFile=adaptorFile;
        this.seedMismatchesTrimmomatic=seedMismatchesTrimmomatic;
        this.palindromeClipThresholdTrimmomatic=palindromeClipThresholdTrimmomatic;
        this.simpleClipThreshold=simpleClipThreshold;
        this.processCutadapt=processCutadapt;
        this.processTrimmomatic=processTrimmomatic;
        this.fastaFileLeftOutlier= new BufferedWriter(new FileWriter(pathFastaFileLeftOutlier));
        this.fastaFileRightOutlier= new BufferedWriter(new FileWriter(pathFastaFileRightOutlier));

    }

    /**
     * Method of the class TrimModes to trim sequence to create sequences files for cutadapt with the side-windows method.
     * @return a HasMap with values update for trimming
     * @throws IOException
     * @throws InterruptedException
     */
    public HashMap trimOutlierWithSWMethod() throws IOException, InterruptedException {

        UtilsTrimming utilsTrimming = new UtilsTrimming();

        Pattern pattern1 = Pattern.compile("(([0-9]*[A-Z]).*)");
        Pattern pattern2 = Pattern.compile("([0-9]*)(.)");
        int longueurSequenceCigar=0;
        BufferedWriter BufferedWriterFastqOutput = new BufferedWriter(new FileWriter(this.nameOutputFastq));
        Utils utils = new Utils();

        for (String ID : this.fastqHash.keySet()) {

            String sequenceCigarBinaire="";
            String[] tabValue = this.fastqHash.get(ID);
            String sequence = tabValue[0];
            String cigar=tabValue[2];
            String score=tabValue[1];


            if(!cigar.equals("*")){

                while(cigar.length()!=0){

                    Matcher matcher1 = pattern1.matcher(cigar);
                    boolean b = matcher1.matches();
                    // si recherche fructueuse
                    if(b) {

                        // pour chaque groupe
                        String oneCigare=matcher1.group(2);
                        cigar=cigar.substring(matcher1.group(2).length(),cigar.length());
                        Matcher matcher2 = pattern2.matcher(oneCigare);
                        boolean b2 = matcher2.matches();
                        // si recherche fructueuse
                        if(b2){

                            longueurSequenceCigar=Integer.parseInt(matcher2.group(1))-1;
                            for(int i=0;i<=longueurSequenceCigar;i++){

                                if(matcher2.group(2).equals("M")){
                                    sequenceCigarBinaire=sequenceCigarBinaire+1;
                                    continue;
                                }
                                if(matcher2.group(2).equals("N") || matcher2.group(2).equals("D")){
                                    continue;
                                }
                                if(!matcher2.group(2).equals("N") || !matcher2.group(2).equals("M")){
                                    sequenceCigarBinaire=sequenceCigarBinaire+0;
                                }
                            }
                        }
                    }
                }

                int lengthOutlierBegin=sideWindowsLeft(sequenceCigarBinaire);
                int lengthOutlierEnd=sideWindowsRight(sequenceCigarBinaire);
                tabValue[3]=""+lengthOutlierBegin;
                tabValue[4]=""+lengthOutlierEnd;

                System.out.println(lengthOutlierBegin+"    "+lengthOutlierEnd+"     "+sequence.length());

                if(this.processCutadapt){
                    utilsTrimming.writeOutliers(lengthOutlierBegin,lengthOutlierEnd,sequence,ID, utils,this.fastaFileLeftOutlier, this.fastaFileRightOutlier);
                }

                if(this.processTrimmomatic){

                    String leftOutlierSequence= utilsTrimming.getOutlierLeftSequence(lengthOutlierBegin,sequence);
                    String rightOutlierSequence= utilsTrimming.getOutlierRightSequence(lengthOutlierEnd,sequence);
                    String leftOutlierScore= utilsTrimming.getOutlierLeftScore(lengthOutlierBegin,score);
                    String rightOutlierScore= utilsTrimming.getOutlierRightScore(lengthOutlierEnd,score);

                    String mainSequence = sequence.substring(lengthOutlierBegin,sequence.length()-lengthOutlierEnd);
                    String mainScore = score.substring(lengthOutlierBegin,score.length()-lengthOutlierEnd);

                    TrimWithTrimmomatic trimmingTrimmomatic= new TrimWithTrimmomatic(this.adaptorFile, this.seedMismatchesTrimmomatic, this.palindromeClipThresholdTrimmomatic, this.simpleClipThreshold);

                    String rigthTrimSequence=trimmingTrimmomatic.TrimmomaticTrim(rightOutlierSequence,rightOutlierScore);
                    String leftTrimSequence=trimmingTrimmomatic.TrimmomaticTrim(utils.reverse(leftOutlierSequence),utils.reverse(leftOutlierScore));

                    String rigthTrimScore=score.substring(mainScore.length(),mainScore.length()+rightOutlierSequence.length());
                    String leftTrimScore=score.substring(leftOutlierScore.length()-leftTrimSequence.length(),leftOutlierScore.length());

                    String sequenceTranscript=leftTrimSequence+mainSequence+rigthTrimSequence;
                    String scoreTranscript = leftTrimScore+mainScore+rigthTrimScore;

                    utils.writeFastq( ID, sequenceTranscript, scoreTranscript, BufferedWriterFastqOutput);

                }
            }
        }
        return this.fastqHash;
    }

    //
    // get index by side window
    //

    /**
     * Method of the class TrimModes to the length of the left outlier with a binair CIGAR sequence.
     * @param sequenceCigarBinaire, a binaire CIGAR sequence
     * @return int, the length of the outlier
     */
    private int sideWindowsLeft(String sequenceCigarBinaire){

        if(sequenceCigarBinaire==null){
            return 0;
        }
        String windows="";
        int length=0;
        for(int i = 0;i<=sequenceCigarBinaire.length();i++){
            if(i==(sequenceCigarBinaire.length()-this.lengthWindowsSW-2)){
                break;
            }
            windows=sequenceCigarBinaire.substring(i,i+this.lengthWindowsSW);
            if(sumWindowCIGAR(windows)>=this.thresholdSW){
                length=i+this.lengthWindowsSW;
                break;
            }
        }
        return length;
    }

    /**
     * Method of the class TrimModes to the length of the right outlier with a binair CIGAR sequence.
     * @param sequenceCigarBinaire, a binaire CIGAR sequence
     * @return int, the length of the outlier
     */
    private int sideWindowsRight(String sequenceCigarBinaire){

        if(sequenceCigarBinaire==null){
            return 0;
        }
        String windows="";
        int length=0;
        for(int i = sequenceCigarBinaire.length();i>=0;i--){
            if(i==this.lengthWindowsSW){
                break;
            }

            windows=sequenceCigarBinaire.substring(i-this.lengthWindowsSW,i);
            if(sumWindowCIGAR(windows)>=this.thresholdSW){
                if(i>=sequenceCigarBinaire.length()-this.lengthWindowsSW){
                    length=sequenceCigarBinaire.length()-i;
                }else{
                    length=sequenceCigarBinaire.length()-i-this.lengthWindowsSW;
                }
                break;
            }
        }
        return length;
    }


    /**
     * Method of the class TrimModes to sumWindowCIGAR a sequence CIGAR encode in 1/0.
     * @param windows, a String of CIGAR sequence encode in 1/0
     * @return double, the sumWindowCIGAR of the CIGAR window
     */
    private double sumWindowCIGAR(String windows){
        String[] arrayWindows = windows.split("");
        int sum=0;
        for(int i = 0;i<=arrayWindows.length-1;i++){
            if(arrayWindows[i].equals("1")){
                sum+=1;
            }
        }
        return (double)sum/windows.length();
    }
}
