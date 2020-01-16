folders= 3
list=newArray(folders)
print(list.length)



for (i=0;i<folders;i++) {
	dir = getDirectory("Choose a Directory ");
	list[i] = dir;
	print(list[i]);
}


Array.print(list)

for(i=0; i<folders; i++){
	run("Image Sequence...", "open=" + list[i]+"shifted_img" + " sort" + " use");
	saveAs("Raw Data", list[i]+"shifted_img.raw");
	print(list[i] + " finished");
	close();
}
