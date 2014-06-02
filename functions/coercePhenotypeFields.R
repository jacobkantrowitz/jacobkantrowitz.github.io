## NEED TO COMMENT ON THIS FUNCTION
coercePhenotypeFields <- function(exprData, type)
{
	# COERCE DATA TO APPROPRIATE TYPES
	lbls <- varLabels(exprData)
	exprData2 <- exprData
	
	# for each phenotype fields in exprData
	for(i in 1:traitNumber(exprData))
	{
		# if the label being checked is the same label in the table at this index
		if(lbls[i]==type[i,2])
		{
			# if the class of the label being checked is different than what it should be
			if(class(pData(exprData)[,i])!=type[i,3])
			{
				# change the class of the label being checked to what it should be
				
				if(type[i,3]=="factor")
				{
					cat("Changing trait ", i, "-", as.character(type[i,2]), " to factor\n", sep="")
					pData(exprData2)[,i] <- as.factor(pData(exprData2)[,i])
				}
				else if(type[i,3]=="numeric")
				{
					cat("Changing trait ", i, "-", as.character(type[i,2]), " to numeric\n", sep="")
					pData(exprData2)[,i] <- as.numeric(pData(exprData2)[,i])
				}
				else if(type[i,3]=="integer")
				{
					cat("Changing trait ", i, "-", as.character(type[i,2]), " to integer\n", sep="")
					pData(exprData2)[,i] <- as.integer(pData(exprData2)[,i])
				}
				else if(type[i,3]=="character")
				{
					cat("Changing trait ", i, "-", as.character(type[i,2]), " to character\n", sep="")
					pData(exprData2)[,i] <- as.character(pData(exprData2)[,i])
				}
				else
				{
					cat("ERROR: no class type name: ", type[i,3], "\n", sep="")
				}
				
			}

		}
	
	}
	return(exprData2)
}	