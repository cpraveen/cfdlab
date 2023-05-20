#!/bin/bash

USERNAME=aadi
FULLNAME="Aadi Bhure"
PASSWORD="abhure123"
SECONDARY_GROUPS="staff" 

# ====

if [[ $UID -ne 0 ]]; then echo "Please run $0 as root." && exit 1; fi

# Find out the next available user ID
MAXID=$(dscl . -list /Users UniqueID | awk '{print $2}' | sort -ug | tail -1)
USERID=$((MAXID+1))

# Create the user account
dscl . -create /Users/$USERNAME
dscl . -create /Users/$USERNAME UserShell /bin/zsh
dscl . -create /Users/$USERNAME RealName "$FULLNAME"
dscl . -create /Users/$USERNAME UniqueID "$USERID"
dscl . -create /Users/$USERNAME PrimaryGroupID 20
dscl . -create /Users/$USERNAME NFSHomeDirectory /Users/$USERNAME

dscl . -passwd /Users/$USERNAME $PASSWORD


# Add use to any specified groups
for GROUP in $SECONDARY_GROUPS ; do
    dseditgroup -o edit -t user -a $USERNAME $GROUP
done

# Create the home directory
createhomedir -c -u $USERNAME > /dev/null

echo "Created user #$USERID: $USERNAME ($FULLNAME)"
