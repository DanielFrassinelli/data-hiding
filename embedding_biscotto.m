function [watermarked] = embedding_biscotto(img, watermark)

watermarked = spread_spectrum(img, watermark, 0.1);

end
