# Distribution of Zulf's Jung Personality Types


nsample <- 1000000
vf <- matrix(rghyp( nsample, bigfiveDist),ncol=5)
jvf <- t(jungTransform( vf ))
typesf <- getPersonalityType( jvf, level = 2, ghfit )

# Check uniqueness achieved here
print(length(unique(typesf)))