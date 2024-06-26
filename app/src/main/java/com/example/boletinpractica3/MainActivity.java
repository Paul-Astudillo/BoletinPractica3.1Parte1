package com.example.boletinpractica3;

import androidx.annotation.Nullable;
import androidx.appcompat.app.AppCompatActivity;
import androidx.core.app.ActivityCompat;

import android.Manifest;
import android.app.ProgressDialog;
import android.content.Intent;
import android.content.pm.PackageManager;
import android.graphics.Bitmap;
import android.graphics.drawable.BitmapDrawable;
import android.net.Uri;
import android.os.Build;
import android.os.Bundle;
import android.provider.MediaStore;
import android.view.View;
import android.widget.Button;
import android.widget.ImageView;
import android.widget.TextView;
import android.widget.Toast;

import com.example.boletinpractica3.databinding.ActivityMainBinding;

import java.io.IOException;

public class MainActivity extends AppCompatActivity {

    private static final int REQUEST_PERMISSION_CAMERA = 200;
    private static final int REQUEST_IMAGE_CAMERA = 201;
    private static final int REQUEST_IMAGE_GALLERY = 202;

    private ImageView foto;
    private Button boton, botonMomentosHu, botonMomentosZernike, botonSelectImage;
    private TextView resultTextView, resultTextView2;
    private ProgressDialog progressDialog;

    static {
        System.loadLibrary("boletinpractica3");
    }

    private ActivityMainBinding binding;

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);

        binding = ActivityMainBinding.inflate(getLayoutInflater());
        setContentView(binding.getRoot());

        foto = findViewById(R.id.imageView);
        boton = findViewById(R.id.button);

        boton.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View v) {
                if (Build.VERSION.SDK_INT >= Build.VERSION_CODES.M) {
                    if (ActivityCompat.checkSelfPermission(MainActivity.this, Manifest.permission.CAMERA) == PackageManager.PERMISSION_GRANTED) {
                        encenderCamara();
                    } else {
                        ActivityCompat.requestPermissions(MainActivity.this, new String[]{Manifest.permission.CAMERA}, REQUEST_PERMISSION_CAMERA);
                    }
                } else {
                    encenderCamara();
                }
            }
        });

        resultTextView = findViewById(R.id.resultTextView);
        botonMomentosHu = findViewById(R.id.button_momentos_hu);
        botonMomentosHu.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View v) {
                calcularMomentosHu();
            }
        });

        resultTextView2 = findViewById(R.id.resultTextView2);
        botonMomentosZernike = findViewById(R.id.button_momentos_zernike);
        botonMomentosZernike.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View v) {
                calcularMomentosZernike();
            }
        });

        botonSelectImage = findViewById(R.id.button_select_image);
        botonSelectImage.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View v) {
                selectImageFromGallery();
            }
        });
    }

    @Override
    public void onRequestPermissionsResult(int requestCode, String[] permissions, int[] grantResults) {
        super.onRequestPermissionsResult(requestCode, permissions, grantResults);
        if (requestCode == REQUEST_PERMISSION_CAMERA) {
            if (grantResults.length > 0 && grantResults[0] == PackageManager.PERMISSION_GRANTED) {
                encenderCamara();
            } else {
                Toast.makeText(this, "Permiso de cámara denegado", Toast.LENGTH_SHORT).show();
            }
        } else if (requestCode == REQUEST_IMAGE_GALLERY) {
            if (grantResults.length > 0 && grantResults[0] == PackageManager.PERMISSION_GRANTED) {
                selectImageFromGallery();
            } else {
                Toast.makeText(this, "Permiso de almacenamiento denegado", Toast.LENGTH_SHORT).show();
            }
        }
    }

    @Override
    protected void onActivityResult(int requestCode, int resultCode, @Nullable Intent data) {
        super.onActivityResult(requestCode, resultCode, data);
        if (requestCode == REQUEST_IMAGE_CAMERA && resultCode == RESULT_OK && data != null) {
            Bitmap imageBitmap = (Bitmap) data.getExtras().get("data");
            foto.setImageBitmap(imageBitmap);
        } else if (requestCode == REQUEST_IMAGE_GALLERY && resultCode == RESULT_OK && data != null) {
            Uri selectedImage = data.getData();
            try {
                Bitmap bitmap = MediaStore.Images.Media.getBitmap(this.getContentResolver(), selectedImage);
                foto.setImageBitmap(bitmap);
            } catch (IOException e) {
                e.printStackTrace();
                Toast.makeText(this, "Error al cargar la imagen", Toast.LENGTH_SHORT).show();
            }
        }
    }

    private void encenderCamara() {
        Intent cameraIntent = new Intent(MediaStore.ACTION_IMAGE_CAPTURE);
        if (cameraIntent.resolveActivity(getPackageManager()) != null) {
            startActivityForResult(cameraIntent, REQUEST_IMAGE_CAMERA);
        }
    }

    private void selectImageFromGallery() {
        Intent intent = new Intent(Intent.ACTION_PICK, MediaStore.Images.Media.EXTERNAL_CONTENT_URI);
        if (intent.resolveActivity(getPackageManager()) != null) {
            startActivityForResult(intent, REQUEST_IMAGE_GALLERY);
        }
    }

    private void calcularMomentosHu() {
        if (foto.getDrawable() == null) {
            Toast.makeText(this, "Primero capture una imagen", Toast.LENGTH_SHORT).show();
            return;
        }

        Bitmap tomada = ((BitmapDrawable) foto.getDrawable()).getBitmap();
        if (tomada == null) {
            Toast.makeText(this, "Error: Imagen vacía", Toast.LENGTH_SHORT).show();
            return;
        }

        progressDialog = ProgressDialog.show(this, "Calculando", "Por favor espere...", true);

        new Thread(new Runnable() {
            @Override
            public void run() {
                String result = calcularHuMoments(tomada);
                runOnUiThread(new Runnable() {
                    @Override
                    public void run() {
                        progressDialog.dismiss();
                        resultTextView.setText(result);
                    }
                });
            }
        }).start();
    }

    private void calcularMomentosZernike() {
        if (foto.getDrawable() == null) {
            Toast.makeText(this, "Primero capture una imagen", Toast.LENGTH_SHORT).show();
            return;
        }

        Bitmap tomada = ((BitmapDrawable) foto.getDrawable()).getBitmap();
        if (tomada == null) {
            Toast.makeText(this, "Error: Imagen vacía", Toast.LENGTH_SHORT).show();
            return;
        }

        progressDialog = ProgressDialog.show(this, "Calculando", "Por favor espere...", true);

        new Thread(new Runnable() {
            @Override
            public void run() {
                String result = calcularZernikeMoments(tomada);
                runOnUiThread(new Runnable() {
                    @Override
                    public void run() {
                        progressDialog.dismiss();
                        resultTextView2.setText(result);
                    }
                });
            }
        }).start();
    }

    public native String calcularHuMoments(Bitmap bitmap);
    public native String calcularZernikeMoments(Bitmap bitmap);
}
