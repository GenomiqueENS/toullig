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
 * Created by birer on 27/03/17.
 */
public class TrimModes {

    private int lengthWindowsSW;
    private double thresholdSW;
    private HashMap<String, String[]> fastqHash;
    private File nameOutputFastq;
    private int addIndexOutlier;
    private File adaptorFile;
    private int seedMismatchesTrimmomatic;
    private int palindromeClipThresholdTrimmomatic;
    private int simpleClipThreshold;
    private boolean processCutadapt;
    private boolean processTrimmomatic;
    private BufferedWriter fastaFileLeftOutlier;
    private BufferedWriter fastaFileRightOutlier;

    public TrimModes(int lengthWindowsSW, double thresholdSW, HashMap fastqHash, File nameOutputFastq, int addIndexOutlier,File adaptorFile, int seedMismatchesTrimmomatic, int palindromeClipThresholdTrimmomatic, int simpleClipThreshold, boolean processCutadapt, boolean processTrimmomatic, String pathFastaFileLeftOutlier, String pathFastaFileRightOutlier) throws IOException {

        this.lengthWindowsSW=lengthWindowsSW;
        this.thresholdSW=thresholdSW;
        this.fastqHash=fastqHash;
        this.nameOutputFastq=nameOutputFastq;
        this.addIndexOutlier=addIndexOutlier;
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
     * Method of the class TrimFastq to trim sequence to create sequences files for cutadapt.
     * @throws IOException
     * @throws InterruptedException
     */
    public HashMap trimOutlierWithPMethod() throws IOException, InterruptedException {

        int count1=0;
        int count2=0;
        int count3=0;
        int count4=0;
        int count5=0;
        int count6=0;
        int count7=0;

        BufferedWriter BufferedWriterFastqOutput = new BufferedWriter(new FileWriter(this.nameOutputFastq));

        for (String ID : this.fastqHash.keySet()) {

            int lengthOutlierEnd=0;
            int lengthOutlierBegin=0;
            count1++;
            String[] tabValue = this.fastqHash.get(ID);
            String sequence = tabValue[0];
            String score = tabValue[1];
            String cigar = tabValue[2];
            //trim by CIGAR
            if(!cigar.equals("*")){

                //Repere Soft and Hard clipping
                int beginIndexExtremite=cigar.indexOf("S");
                if(beginIndexExtremite>cigar.indexOf("H") && tabValue[2].indexOf("H")>=0){
                    beginIndexExtremite=cigar.indexOf("H");
                }

                if(beginIndexExtremite!=-1 && beginIndexExtremite+1!=cigar.length() && beginIndexExtremite<=6){

                    lengthOutlierBegin = Integer.parseInt(cigar.substring(0,beginIndexExtremite))+this.addIndexOutlier;
                    tabValue[3]=""+lengthOutlierBegin;

                    count2++;
                }else{
                    //System.out.println("CIGAR problem deb: "+tabValue[2]);
                    lengthOutlierBegin=0;
                }

                //Repere Soft and Hard clipping
                int endIndexExtremite=cigar.lastIndexOf("S");
                if(endIndexExtremite<cigar.lastIndexOf("H")){
                    endIndexExtremite=cigar.lastIndexOf("H");
                }
                if(endIndexExtremite==-1){
                    endIndexExtremite=0;
                }

                if(endIndexExtremite+1==cigar.length() && cigar.length()-endIndexExtremite<5 ){

                    lengthOutlierEnd = Integer.parseInt(cigar.substring(cigar.lastIndexOf("M")+1,endIndexExtremite))+this.addIndexOutlier;
                    tabValue[4]=""+lengthOutlierEnd;
                    count3++;
                }else{
                    //System.out.println("CIGAR problem end: "+tabValue[2]);
                    lengthOutlierEnd=0;
                }

                if(this.processCutadapt){
                    writeOutlier(lengthOutlierBegin,lengthOutlierEnd,sequence,ID);
                }

                if(this.processTrimmomatic){

                    Utils utils = new Utils();

                    String leftOutlierSequence= getOutlierLeftSequence(lengthOutlierBegin,sequence);
                    String rightOutlierSequence= getOutlierRightSequence(lengthOutlierEnd,sequence);
                    String leftOutlierScore= getOutlierLeftScore(lengthOutlierBegin,score);
                    String rightOutlierScore= getOutlierRightScore(lengthOutlierEnd,score);

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
            }else{
                count4++;
            }
            if(tabValue[5].equals("16")){
                count5++;
            }
            if(tabValue[5].equals("0")){
                count6++;
            }
            if(tabValue[5].equals("4")){
                count7++;
            }
        }
        System.out.println("Nombre total de sequence dans le SAM : "+count1);
        System.out.println("Nombre de sequence CIGAR '*': "+count4);
        System.out.println("Nombre de QFlag '16': "+count5);
        System.out.println("Nombre de QFlag '0': "+count6);
        System.out.println("Nombre de QFlag '4': "+count7);
        System.out.println("Nombre d'Outlier Debut trouvé: "+count2);
        System.out.println("Nombre d'Outlier Fin trouvé: "+count3);
        if(this.processTrimmomatic) {
            BufferedWriterFastqOutput.close();
        }
        return this.fastqHash;
    }

    /**
     * Method of the class TrimFastq to trim sequence to create sequences files for cutadapt with the side-windows method.
     * @throws IOException
     * @throws InterruptedException
     */
    public HashMap trimOutlierWithSWMethod() throws IOException, InterruptedException {

        Pattern pattern1 = Pattern.compile("(([0-9]*[A-Z]).*)");
        Pattern pattern2 = Pattern.compile("([0-9]*)(.)");
        int longueurSequenceCigar=0;
        BufferedWriter BufferedWriterFastqOutput = new BufferedWriter(new FileWriter(this.nameOutputFastq));

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
                    writeOutlier(lengthOutlierBegin,lengthOutlierEnd,sequence,ID);
                }

                if(this.processTrimmomatic){

                    Utils utils = new Utils();

                    String leftOutlierSequence= getOutlierLeftSequence(lengthOutlierBegin,sequence);
                    String rightOutlierSequence= getOutlierRightSequence(lengthOutlierEnd,sequence);
                    String leftOutlierScore= getOutlierLeftScore(lengthOutlierBegin,score);
                    String rightOutlierScore= getOutlierRightScore(lengthOutlierEnd,score);

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
    // get Outlier
    //

    /**
     * Method of the class TrimFastq to write the outliers (3' and 5') in fasta files.
     * @param lengthOutlierBegin, the length of the outlier right
     * @param lengthOutlierEnd, the length of the outlier left
     * @param sequence, the sequence of the read
     * @param ID, the ID of the read
     * @throws IOException
     */
    private void writeOutlier(int lengthOutlierBegin, int lengthOutlierEnd, String sequence, String ID) throws IOException {

        String leftOutlierSequence= getOutlierLeftSequence(lengthOutlierBegin,sequence);
        String rightOutlierSequence= getOutlierRightSequence(lengthOutlierEnd,sequence);
        Utils utils = new Utils();
        utils.writeFasta(leftOutlierSequence,ID,this.fastaFileLeftOutlier);
        utils.writeFasta(rightOutlierSequence,ID,this.fastaFileRightOutlier);
    }

    /**
     * Method of the class TrimFastq to obtain the sequence of the left outlier.
     * @param lengthOutlierBegin, the length of the outlier
     * @param sequence, the sequence of the read
     * @return the sequence of the left outlier
     */
    private String getOutlierLeftSequence(int lengthOutlierBegin, String sequence){
        return sequence.substring(0,lengthOutlierBegin);
    }

    /**
     * Method of the class TrimFastq to obtain the sequence of the right outlier.
     * @param lengthOutlierEnd, the length of the outlier
     * @param sequence, the sequence of the read
     * @return the sequence of the right outlier
     */
    private String getOutlierRightSequence(int lengthOutlierEnd, String sequence){

        return sequence.substring(sequence.length()-lengthOutlierEnd,sequence.length());
    }

    /**
     * Method of the class TrimFastq to obtain the score of the left outlier.
     * @param lengthOutlierBegin, the length of the outlier
     * @param score, the score of the read
     * @return
     */
    private String getOutlierLeftScore(int lengthOutlierBegin, String score){
        return score.substring(0,lengthOutlierBegin);
    }

    /**
     * Method of the class TrimFastq to obtain the score of the right outlier.
     * @param lengthOutlierEnd, the length of the outlier
     * @param score, the score of the read
     * @return
     */
    private String getOutlierRightScore(int lengthOutlierEnd, String score){
        return score.substring(score.length()-lengthOutlierEnd,score.length());
    }

    //
    // get index by side window
    //

    /**
     * Method of the class TrimFastq to the length of the left outlier with a binair CIGAR sequence.
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
     * Method of the class TrimFastq to the length of the right outlier with a binair CIGAR sequence.
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
     * Method of the class TrimFastq to sumWindowCIGAR a sequence CIGAR encode in 1/0.
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
