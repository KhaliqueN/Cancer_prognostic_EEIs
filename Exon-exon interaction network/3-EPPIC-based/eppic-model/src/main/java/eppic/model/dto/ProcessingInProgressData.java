package eppic.model.dto;


/**
 * This class is used to transfer information necessary to display the status of submitted job.
 * @author srebniak_a
 *
 */
public class ProcessingInProgressData implements ProcessingData 
{
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	/**
	 * Identitifier of the job.
	 */
	private String jobId;
	
	/**
	 * Current status of the job.
	 */
	private String status;
	
	/**
	 * Output of the job processing.
	 */
	private String log;
	
	/**
	 * Pdb code or name of the file submitted.
	 */
	private String inputName;
	
	/**
	 * Type of the inputName - pdb code or file.
	 */
	private int inputType;
	
	/**
	 * Current step.
	 */
	private StepStatus step;

	/**
	 * Retrieves identifier of the job.
	 * @return identifier of the job
	 */
	public String getJobId() {
		return jobId;
	}

	/**
	 * Sets identifier of the job.
	 * @param jobId identifier of the job
	 */
	public void setJobId(String jobId) {
		this.jobId = jobId;
	}

	/**
	 * Retrieves current status of the job.
	 * @return current status of the job
	 */
	public String getStatus() {
		return status;
	}

	/**
	 * Sets current status of the job.
	 * @param status current status of the job
	 */
	public void setStatus(String status) {
		this.status = status;
	}

	/**
	 * Retrieves output of processing.
	 * @return output of processing
	 */
	public String getLog() {
		return log;
	}

	/**
	 * Sets output of processing
	 * @param log output of processing
	 */
	public void setLog(String log) {
		this.log = log;
	}

	/**
	 * Sets pdb code or name of submitted file.
	 * @param inputName pdb code or name of submitted file
	 */
	public void setInputName(String inputName) {
		this.inputName = inputName;
	}

	/**
	 * Retrieves pdb code or name of submitted file
	 * @return pdb code or name of submitted file
	 */
	public String getInputName() {
		return inputName;
	}

	/**
	 * Sets current step.
	 * @param step current step
	 */
	public void setStep(StepStatus step) {
		this.step = step;
	}

	/**
	 * Retrieves current step.
	 * @return current step
	 */
	public StepStatus getStep() {
		return step;
	}

	/**
	 * Sets type of the inputName
	 * @param inputType type of the inputName
	 */
	public void setInputType(int inputType) {
		this.inputType = inputType;
	}

	/**
	 * Retrieves type of the inputName.
	 * @return type of the inputName
	 */
	public int getInputType() {
		return inputType;
	}

}
