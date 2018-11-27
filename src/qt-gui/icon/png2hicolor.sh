mkdir -p hicolor

for i in 128 64 48 32 22 16; do
  mkdir -p hicolor/${i}x${i}/apps
  convert icon.png -resize ${i}x${i} hicolor/${i}x${i}/apps/maracluster.png
done
