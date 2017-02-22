package fr.ens.biologie.genomique.toullig;

import java.util.*;

public class LocalReportertoLog {

    private List<String> listWriteSequenceLog = new ArrayList<String>();
    private List<String> listWorkflowStatusLog = new ArrayList<String>();

    private LocalReporter localReporter = new LocalReporter();
    private Fast5toFastq f5;

    public LocalReportertoLog(Fast5toFastq runConversionFastq){
        this.localReporter=runConversionFastq.getLocalReporter();
        this.f5 = runConversionFastq;
        getAllListForLog(localReporter);
    }

    //
    //
    // Create Log
    //
    //

    /**
     * This method of the class Fast5toFastq get a list of log on the number of
     * files.
     * @return a list of string
     */
    public List<String> getListLog() {
        List<String> listLog = new ArrayList<String>();
        listLog.add("Number of fast5 files read " + this.f5.getNumberFast5Files(this.localReporter));
        listLog.add("Number of corrupt files " + this.f5.getNumberCorruptFast5Files(this.localReporter));
        listLog.add("Number of calibrate strand files read "
                + this.f5.getNumberCalibrateStrandFast5Files(this.localReporter));
        listLog.add("Number of unclassified files read "
                + this.f5.getNumberUnclassifiedFast5Files(this.localReporter));
        listLog.add("Number of fail files read " + this.f5.getNumberFailFast5Files(this.localReporter));
        listLog.add("Number of pass files read " + this.f5.getNumberPassFast5Files(this.localReporter));

        listLog.add(
                "Number of pass barcode file " + this.f5.getNumberBarcodeFast5Files(this.localReporter) + "\n");
        for (String element : this.listWriteSequenceLog) {
            listLog.add(element);
        }
        return listLog;
    }

    /**
     * This method of the class Fast5toFastq get a list of log on the status of
     * the workflows.
     * @return a list of string
     */
    public List<String> getListLogStatusWorkflow() {
        return this.listWorkflowStatusLog;
    }

    /**
     * This method of the class Fast5toFastq get All Log for the conversion of
     * fast5 to fastq file.
     * @param localReporter a LocalReporter
     */
    private void getAllListForLog(LocalReporter localReporter) {

        //
        // Count of write sequence in fastq file
        //
        for(String element : localReporter.getCounterNames("numberSequenceWrite")){
            long numberSequence= localReporter.getCounterValue("numberSequenceWrite", element);
            if(numberSequence==-1){
                        numberSequence=0;
            }
            String[] partSplit = element.split("_");
            this.listWriteSequenceLog.add("In the file "
                    + partSplit[0] + " template the number of total sequence write "
                    + numberSequence);
        }
        List<String> groupName=new ArrayList<>();
        for(String element : localReporter.getCounterGroups()){
            groupName.add(element);
        }

        Collections.sort(groupName);
        String oldStatus="";
        for(String element : groupName){

            String[] part=element.split("_");
            if(element.contains("_")){
                String status = part[0];
                if(!oldStatus.equals(status)){
                    this.listWorkflowStatusLog.add("\n###################################\n"
                            + status + " :\n###################################");
                }

                this.listWorkflowStatusLog.add("\n"+part[1]+" :\n");
                stackHashWorkflow(localReporter, element, part[0]);

                oldStatus=status;
            }

        }

    }


    /**
     * This method of the class Fast5toFastq stack the list of workflow status.
     * @param localReporter, name of object LocalReporter
     * @param group, name of group
     * @param status, the name of folder fast5 status
     */
    private void stackHashWorkflow(LocalReporter localReporter, String group,
                                   String status) {
        Set<String> keys = localReporter.getCounterNames(group);

        if (keys.size() >= 1) {
            for (String key : keys) {

                if(status.equals("fail") && key.contains("Calibration strand detected.") ){
                    localReporter.incrCounter("numberFiles",
                            "numberCalibrateStrandFast5Files", 1);
                }
                this.listWorkflowStatusLog.add("The status \""
                        + key + "\" is present " + localReporter.getCounterValue(group, key)
                        + " times in the folder " + status);
            }
        } else {
            this.listWorkflowStatusLog.add("No status for the "
                    + group + " Workflow in the folder " + status);
        }
    }
}
