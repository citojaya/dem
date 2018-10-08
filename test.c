#include <stdio.h>

int main(){
    FILE *fp;
    char buff[255];

    fp = fopen("text.txt","r");
    // fscanf(fp,"%s",buff);
    // printf("%s\n",buff);

    // fgets(buff,10,(FILE*)fp);
    // printf("2: %s\n",buff);

    // fgets(buff,255,(FILE*)fp);
    // printf("3: %s\n",buff);

    // fgets(buff,255,(FILE*)fp);
    // printf("4: %s\n",buff);
    char str1[10],str2[10];
    //,str2[10],str3[10],str4[10];
    fscanf(fp,"%s %s", str1, str2);
    printf("%s %s",str1, str2);
    
    fclose(fp);

}