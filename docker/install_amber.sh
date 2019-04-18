cd amber18
export AMBERHOME=`pwd`
./configure gnu
     # We recommend you say "yes" when asked to apply updates
source amber.sh # Use amber.csh if you use tcsh or csh
make install
echo "source $AMBERHOME/amber.sh" >> ~/.bashrc  # Add Amber to your environment
