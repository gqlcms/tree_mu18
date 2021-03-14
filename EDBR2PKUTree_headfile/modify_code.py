import re

# a = re.compile('fChain->SetBranchAddress("^[\w\s\S&]$"^[\w\s\S&,]$);')
a = re.compile('SetBranchAddress')
num = 10
count = 0
for i in open("branch_need_toSet.Cpp","r").readlines():
    if count<num:
        if 'SetBranchAddress' in i:
            new_text = a.sub(r'hhh',i)
            print new_text
            count += 1

text='hhhh 12/18/19'
print(re.sub(r'(\d+)/(\d+)/(\d+)',r'\3-\1-\2',text))
print(text)
    